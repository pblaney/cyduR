#!/usr/bin/env Rscript

# Execute the cyduR workflow

#########################
#####   Libraries   #####

suppressPackageStartupMessages({
  library(optparse)
  library(cli)
  library(clisymbols)
})

#########################
#####   Functions   #####

# -----------------------------------------------------------
# Error Handler
# -----------------------------------------------------------

handle_error <- function(e, context = "Script error") {
  cli_div(theme = list(
    ".error-box" = list(
      "padding-left" = 2,
      "border-left" = "2px solid #FF5555"
    )
  ))
  cli_div(class = "error-box")
  cli_alert_danger("{context}")
  cli_text("{.strong Message:} {e$message}")
  cli_end()
  quit(save = "no", status = 1)
}

handle_warning <- function(w) {
  cli_alert_warning("{w$message}")
  invokeRestart("muffleWarning")
}

safe_run <- function(expr, context = "Executing task") {
  tryCatch(
    eval(substitute(expr)),
    error = function(e) handle_error(e, context),
    warning = function(w) handle_warning(w)
  )
}

# Read a FASTA file
read_fasta <- function(fasta_file) {
  
  cli_alert_info("Parsing FASTA sequences")
  
  # Read all lines
  lines <- readLines(fasta_file, warn = FALSE)
  
  names_list <- c()
  data <- list()
  
  current_id <- NULL
  current_seq <- NULL
  idx <- 0
  
  for (line in lines) {
    
    # Header line
    if (startsWith(line, ">")) {
      
      # If we were recording a previous sequence, store it before starting a new one
      if (!is.null(current_id)) {
        idx <- idx + 1
        names_list[idx] <- current_id
        data[[idx]] <- current_seq
        names(data)[idx] <- current_id
      }
      
      # Start new sequence
      current_id <- sub("^>", "", line)
      current_seq <- ""
      
    } else {
      # Append sequence line
      current_seq <- paste0(current_seq, line)
    }
  }
  
  # Store the final sequence
  if (!is.null(current_id)) {
    idx <- idx + 1
    names_list[idx] <- current_id
    data[[idx]] <- current_seq
    names(data)[idx] <- current_id
  }
  
  return(list(names_list = names_list, data = data))
}

# Helper function to open and append FASTA files
append_fasta_entry <- function(fasta_file, entry) {
  
  con <- file(fasta_file, open = "a")
  on.exit(close(con), add = TRUE)
  
  #--------------------------
  # Validate file
  #--------------------------
  if (!file.exists(fasta_file)) {
    stop("The FASTA file does not exist: ", fasta_file)
  }
  
  #--------------------------
  # Validate entry
  #--------------------------
  if (!is.character(entry) || length(entry) == 0) {
    stop("`entry` must be a non-empty character string or vector.")
  }
  
  # Collapse vector input into a single string, then re-split into lines
  entry <- paste(entry, collapse = "\n")
  entry_lines <- unlist(strsplit(entry, "\n"))
  
  # Header must begin with ">"
  header <- entry_lines[1]
  if (!grepl("^>", header)) {
    stop("FASTA entry must begin with a header line starting with '>'.")
  }
  
  # Ensure a sequence exists
  if (length(entry_lines) < 2) {
    stop("The FASTA entry must contain at least one sequence line.")
  }
  
  # Remove whitespace inside sequence lines
  seq_lines <- entry_lines[-1]
  seq_lines <- gsub("\\s+", "", seq_lines)
  
  # Final lines to write
  output_lines <- c(header, seq_lines)
  
  #--------------------------
  # Append to FASTA file
  #--------------------------
  writeLines(output_lines, con = con)
  
  invisible(TRUE)
}

# This function creates a sequence with GC-content close to the GC-content of the input sequence,
# but the counts of nucleotides may differ from the input sequence.
gc3 <- function(seq) {
  
  # Initialize the GC and AT content
  gc <- 0
  at <- 0
  
  # First, calculate the A+T and G+C pairs of the input sequence in third codon position
  # Iterate through 3rd codon positions in sequence
  for (num in seq(3, nchar(seq), by = 3)) {
    base <- substr(seq, num, num)
    if (base %in% c("A", "T")) {
      at <- at + 1
    } else if (base %in% c("G", "C")) {
      gc <- gc + 1
    }
  }
  
  # Normalize GC and AT counts
  codon_count <- nchar(seq) / 3
  gc <- gc / codon_count
  at <- at / codon_count
  
  # Two part sequence shuffling
  #    1. Flagging for type of randomization based on third codon position
  #    2. Execute the shuffling to build final sequence
  
  # Build seq1 (list of codon bases and flags for randomization)
  seq1 <- c()
  
  # list "seq1" will contain the first two nt of codon, third codon position will
  # contain flags for subsequent randomization. Flags ('_Y_','_R_','_H_',or '_N_')
  # correspond to IUPAC single-letter code, Y-Pyrimindine(C or T), R-Purine(A or G),
  # H-Not G(A or C or T), N-any.
  for (num in seq(3, nchar(seq), by = 3)) {
    third_base <- substr(seq, num, num)
    codon <- substr(seq, num - 2, num - 1)
    
    seq1 <- c(seq1, codon) # Add codon to seq1
    
    if (third_base %in% c("T", "C") && codon %in% c("TT", "TA", "TG", "CA", "AA", "AG", "GA")) {
      seq1 <- c(seq1, "_Y_") # Pyrimidine (C or T)
      
    } else if (third_base %in% c("A", "G") && codon %in% c("TT", "CA", "AA", "AG", "GA")) {
      seq1 <- c(seq1, "_R_") # Purine (A or G)
      
    } else if (substr(seq, num - 2, num + 1) %in% c("ATT", "ATC", "ATA")) {
      seq1 <- c(seq1, "_H_") # Not G (A, C, or T)
      
    } else if (third_base %in% c("A", "G", "T", "C") && codon %in% c("TC", "CT", "CC", "CG", "AC", "GT", "GC", "GG")) {
      seq1 <- c(seq1, "_N_") # Any (A, T, C, or G)
      
    } else {
      seq1 <- c(seq1, third_base) # Default: keep base
    }
  }
  
  # Build seq2 (derived sequence based on GC-content randomization)
  # "seq2" will contain the derived sequence, appropriate nucleotide is chosen for
  # flags in "seq1", according to GC-content
  seq2 <- ""
  for (i in seq1) {
    x <- runif(1) # Random number generation (between 0 and 1)
    
    # Pyrimidine (C or T) 
    if (i == "_Y_") {
      # GC content test
      seq2 <- paste0(seq2, ifelse(x <= gc, "C", ifelse(x <= gc + at, "T", sample(c("T", "C"), 1))))
      
      # Purine (A or G)
    } else if (i == "_R_") {
      # GC content test
      seq2 <- paste0(seq2, ifelse(x <= gc, "G", ifelse(x <= gc + at, "A", sample(c("A", "G"), 1))))
      
    } else if (i == "_H_") {
      # GC content test
      seq2 <- paste0(seq2, ifelse(x <= gc, "C", ifelse(x <= gc + at, sample(c("A", "T"), 1), sample(c("A", "T", "C"), 1))))
      
    } else if (i == "_N_") {
      # GC content test
      seq2 <- paste0(seq2, ifelse(x <= gc, sample(c("G", "C"), 1), ifelse(x <= gc + at, sample(c("A", "T"), 1), sample(c("A", "G", "T", "C"), 1))))
      
    } else {
      seq2 <- paste0(seq2, i) # Add non-flag bases directly
    }
  }
  return(seq2)
}

# This function creates scrambled sequence with the numbers of each nucleotide
# identical to the input sequence.
n3 <- function(seq) {
  # Helper function to randomly permute elements of a vector
  shuffle <- function(vec) {
    return(sample(vec))
  }
  
  # Initialize lists for flags
  Y <- c()
  R <- c()
  H <- c()
  N <- c()
  
  # -----------------------------------------------------------
  # FLAG `_Y_` (Shuffling where T -> C or C -> T maintains amino acid sequence)
  # -----------------------------------------------------------
  seq1 <- c()
  for (num in seq(3, nchar(seq), by = 3)) {
    codon <- substr(seq, num - 2, num - 1)  # Get codon
    third_base <- substr(seq, num, num) # Get third nucleotide
    # Check if shuffling is possible and add `_Y_` flag
    if (third_base %in% c("T", "C") && codon %in% c("TT", "TC", "TA", "TG", "CT", "CC", "CA", "CG", "AT", "AC", "AA", "AG", "GT", "GC", "GA", "GG")) {
      Y <- c(Y, third_base)  # Add third nucleotide to `Y`
      seq1 <- c(seq1, codon, "_Y_") # Append the codon and `_Y_` flag
    } else {
      seq1 <- c(seq1, substr(seq, num - 2, num))  # Add unflagged codon directly
    }
  }
  # now "seq1" contains flag '_Y_' in the third position of all codons, where C->T
  # or T->C shuffling preserves the amino acid sequence (i.e. PHE, SER etc.).
  # #C and T from the original sequence in this case would be extracted into list "Y"
  
  # Shuffle `Y` and replace `_Y_` flag with shuffled values
  Y <- shuffle(Y) # Shuffle the `Y` list
  seq2 <- ""
  for (i in seq_along(seq1)) {
    if (seq1[i] == "_Y_") {
      seq2 <- paste0(seq2, Y[1]) # Replace `_Y_` with a shuffled value from `Y`
      Y <- Y[-1]                 # Remove the used value from `Y`
    } else {
      seq2 <- paste0(seq2, seq1[i]) # Append unmodified codon/nucleotide
    }
  }
  seq <- seq2 # Update the sequence
  
  
  # -----------------------------------------------------------
  # FLAG `_R_` (Shuffling where A -> G or G -> A maintains protein sequence)
  # -----------------------------------------------------------
  seq1 <- c()
  for (num in seq(3, nchar(seq), by = 3)) {
    codon <- substr(seq, num - 2, num -1)  # Get codon
    third_base <- substr(seq, num, num) # Get third nucleotide
    # Check if shuffling is possible and add `_R_` flag
    if (third_base %in% c("A", "G") && codon %in% c("TT", "TC", "CT", "CC", "CA", "CG", "AC", "AA", "AG", "GT", "GC", "GA", "GG")) {
      R <- c(R, third_base)  # Add third nucleotide to `R`
      seq1 <- c(seq1, codon, "_R_") # Append the codon and `_R_` flag
    } else {
      seq1 <- c(seq1, substr(seq, num - 2, num))  # Add unflagged codon directly
    }
  }
  
  # Shuffle `R` and replace `_R_` flag with shuffled values
  R <- shuffle(R)
  seq2 <- ""
  for (i in seq_along(seq1)) {
    if (seq1[i] == "_R_") {
      seq2 <- paste0(seq2, R[1]) # Replace `_R_` with a shuffled value from `R`
      R <- R[-1]                 # Remove the used value from `R`
    } else {
      seq2 <- paste0(seq2, seq1[i]) # Append unmodified codon/nucleotide
    }
  }
  seq <- seq2 # Update the sequence
  
  # -----------------------------------------------------------
  # FLAG `_H_` (Shuffle A, C, T where possible)
  # -----------------------------------------------------------
  seq1 <- c()
  for (num in seq(3, nchar(seq), by = 3)) {
    codon <- substr(seq, num - 2, num - 1)  # Get codon
    third_base <- substr(seq, num, num) # Get third nucleotide
    # Check if shuffling is possible and add `_H_` flag
    if (third_base %in% c("A", "C", "T") && codon %in% c("TC", "CT", "CC", "CG", "AT", "AC", "GT", "GC", "GG")) {
      H <- c(H, third_base)  # Add third nucleotide to `H`
      seq1 <- c(seq1, codon, "_H_") # Append the codon and `_H_` flag
    } else {
      seq1 <- c(seq1, substr(seq, num - 2, num))  # Add unflagged codon directly
    }
  }
  
  # Shuffle `H` and replace `_H_` flag with shuffled values
  H <- shuffle(H)
  seq2 <- ""
  for (i in seq_along(seq1)) {
    if (seq1[i] == "_H_") {
      seq2 <- paste0(seq2, H[1]) # Replace `_H_` with a shuffled value from `H`
      H <- H[-1]                 # Remove the used value from `H`
    } else {
      seq2 <- paste0(seq2, seq1[i]) # Append unmodified codon/nucleotide
    }
  }
  seq <- seq2 # Update the sequence
  
  # -----------------------------------------------------------
  # FLAG `_N_` (Shuffle all four nucleotides: A, C, T, G)
  # -----------------------------------------------------------
  seq1 <- c()
  for (num in seq(3, nchar(seq), by = 3)) {
    codon <- substr(seq, num - 2, num - 1)  # Get codon
    third_base <- substr(seq, num, num) # Get third nucleotide
    # Check if shuffling is possible and add `_N_` flag
    if (third_base %in% c("A", "C", "T", "G") && codon %in% c("TC", "CT", "CC", "CG", "AC", "GT", "GC", "GG")) {
      N <- c(N, third_base)  # Add third nucleotide to `N`
      seq1 <- c(seq1, codon, "_N_") # Append the codon and `_N_` flag
    } else {
      seq1 <- c(seq1, substr(seq, num - 2, num))  # Add unflagged codon directly
    }
  }
  
  # Shuffle `N` and replace `_N_` flag with shuffled values
  N <- shuffle(N)
  seq2 <- ""
  for (i in seq_along(seq1)) {
    if (seq1[i] == "_N_") {
      seq2 <- paste0(seq2, N[1]) # Replace `_N_` with a shuffled value from `N`
      N <- N[-1]                 # Remove the used value from `N`
    } else {
      seq2 <- paste0(seq2, seq1[i]) # Append unmodified codon/nucleotide
    }
  }
  seq <- seq2 # Update the sequence
  
  return(seq)
}

# This function creates a randomized sequence, with dinucleotide frequences in codon
# position 2-3 close to those of the input sequence (not exact since this only counts
# dinucleotide frequency in second and third positions).
dn23 <- function(seq) {
  # Initialize dinucleotide counts
  aa <- ag <- ac <- at <- ga <- gg <- gc <- gt <- 0
  ca <- cg <- cc <- ct <- ta <- tg <- tc <- tt <- 0
  
  # -----------------------------------------------------------
  # Step 1: Calculate dinucleotide frequencies in codon positions 2-3
  # -----------------------------------------------------------
  for (num in seq(3, nchar(seq), by = 3)) {
    second_base <- substr(seq, num - 1, num - 1)
    third_base <- substr(seq, num, num)
    
    if (second_base == "A") {
      if (third_base == "A") {
        aa <- aa + 1
      } else if (third_base == "G") {
        ag <- ag + 1
      } else if (third_base == "C") {
        ac <- ac + 1
      } else if (third_base == "T") {
        at <- at + 1
      }
    } else if (second_base == "G") {
      if (third_base == "A") {
        ga <- ga + 1
      } else if (third_base == "G") {
        gg <- gg + 1
      } else if (third_base == "C") {
        gc <- gc + 1
      } else if (third_base == "T") {
        gt <- gt + 1
      }
    } else if (second_base == "C") {
      if (third_base == "A") {
        ca <- ca + 1
      } else if (third_base == "G") {
        cg <- cg + 1
      } else if (third_base == "C") {
        cc <- cc + 1
      } else if (third_base == "T") {
        ct <- ct + 1
      }
    } else if (second_base == "T") {
      if (third_base == "A") {
        ta <- ta + 1
      } else if (third_base == "G") {
        tg <- tg + 1
      } else if (third_base == "C") {
        tc <- tc + 1
      } else if (third_base == "T") {
        tt <- tt + 1
      }
    }
  }
  
  # Normalize dinucleotide frequencies
  codon_count <- nchar(seq) / 3
  aa <- aa / codon_count
  ag <- ag / codon_count
  ac <- ac / codon_count
  at <- at / codon_count
  ga <- ga / codon_count
  gg <- gg / codon_count
  gc <- gc / codon_count
  gt <- gt / codon_count
  ca <- ca / codon_count
  cg <- cg / codon_count
  cc <- cc / codon_count
  ct <- ct / codon_count
  ta <- ta / codon_count
  tg <- tg / codon_count
  tc <- tc / codon_count
  tt <- tt / codon_count
  
  # -----------------------------------------------------------
  # Step 2: Substitute codons based on dinucleotide probabilities
  # -----------------------------------------------------------
  seq2 <- ""
  for (num in seq(3, nchar(seq), by = 3)) {
    # Build codon sequence
    seq2 <- paste0(seq2, substr(seq, num - 2, num-1)) # Append 1st and 2nd base
    
    # Substitute according to dinucleotide probabilities
    second_base <- substr(seq, num - 1, num - 1)
    third_base <- substr(seq, num, num)
    
    if (second_base == "A" && !substr(seq, num - 2, num) %in% c("TAA", "TAG")) {
      if (third_base %in% c("T", "C")) {
        space <- at + ac
        AT <- at / space
        AC <- ac / space
        x <- runif(1)
        if (x <= AT) {
          seq2 <- paste0(seq2, "T")
        } else if (x <= AT + AC) {
          seq2 <- paste0(seq2, "C")
        }
      } else if (third_base %in% c("A", "G")) {
        space <- aa + ag
        AA <- aa / space
        AG <- ag / space
        x <- runif(1)
        if (x <= AA) {
          seq2 <- paste0(seq2, "A")
        } else if (x <= AA + AG) {
          seq2 <- paste0(seq2, "G")
        }
      } else {
        seq2 <- paste0(seq2, third_base)
      }
      
    } else if (second_base == "G" && !substr(seq, num - 2, num) %in% c("TGA", "TGG")) {
      if(substr(seq, num - 2, num - 2) %in% c("T", "A") && substr(seq, num, num) %in% c("C", "T")) {
        space <- gt + gc
        GT <- gt / space
        GC <- gc / space
        x <- runif(1)
        if (x <= GA) {
          seq2 <- paste0(seq2, "T")
        } else if (x <= GT + GC) {
          seq2 <- paste0(seq2, "C")
        }
      } else if (substr(seq, num - 2, num) %in% c("AGA", "AGG")) {
        space <- ga + gg
        GA <- ga / space
        GG <- gg / space
        x <- runif(1)
        if (x <= GA) {
          seq2 <- paste0(seq2, "A")
        } else if (x <= GA + GG) {
          seq2 <- paste0(seq2, "G")
        }
      } else if (substr(seq, num - 2, num - 2) %in% c("C", "G")) {
        space <- ga + gg + gc + gt
        GA <- ga / space
        GG <- gg / space
        GC <- gc / space
        GT <- gt / space
        x <- runif(1)
        if (x <= GA) {
          seq2 <- paste0(seq2, "A")
        } else if (x <= GA + GG) {
          seq2 <- paste0(seq2, "G")
        } else if (x <= GA + GG + GC) {
          seq2 <- paste0(seq2, "C")
        } else {
          seq2 <- paste0(seq2, "T")
        }
      } else {
        seq2 <- paste0(seq2, third_base)
      }
      
    } else if (second_base == "C") {
      space <- ca + cg + cc + ct
      CA <- ca / space
      CG <- cg / space
      CC <- cc / space
      CT <- ct / space
      x <- runif(1)
      if (x <= CA) {
        seq2 <- paste0(seq2, "A")
      } else if (x <= CA + CG) {
        seq2 <- paste0(seq2, "G")
      } else if (x <= CA + CG + CC) {
        seq2 <- paste0(seq2, "C")
      } else {
        seq2 <- paste0(seq2, "T")
      }
      
    } else if (second_base == "T") {
      if (substr(seq, num - 2, num) %in% c("TTT", "TTC")) {
        space <- tt + tc
        TT <- tt / space
        TC <- tc / space
        x <- runif(1)
        if (x <= TT) {
          seq2 <- paste0(seq2, "T")
        } else if (TT < x && x <= TT+TC) {
          seq2 <- paste0(seq2, "C")
        }
      } else if (substr(seq, num - 2, num) %in% c("TTA", "TTG")) {
        space <- ta + tg
        TA <- ta / space
        TG <- tg / space
        x <- runif(1)
        if (x <= TA) {
          seq2 <- paste0(seq2, "A")
        } else if (TA < x && x <= TA+TG) {
          seq2 <- paste0(seq2, "C")
        }
      } else if (substr(seq, num - 2, num) %in% c("ATT", "ATC", "ATA")) {
        space <- tt + tc + ta
        TT <- tt / space
        TC <- tc / space
        TA <- ta / space
        x <- runif(1)
        
        if (x <= TA) {
          seq2 <- paste0(seq2, "A")
        } else if (x > TA && x <= TA + TC) {
          seq2 <- paste0(seq2, "C")
        } else if (x > TA + TC && x <= TA + TC + TT) {
          seq2 <- paste0(seq2, "T")
        }
      } else if (substr(seq, num - 2, num - 2) %in% c("C", "G")) {
        space <- ta + tg + tc + tt
        TA <- ta / space
        TG <- tg / space
        TC <- tc / space
        TT <- tt / space
        x <- runif(1)
        if (x <= TA) {
          seq2 <- paste0(seq2, "A")
        } else if (x > TA && x <= TA + TG) {
          seq2 <- paste0(seq2, "G")
        } else if (x > TA + TG && x <= TA + TG + TC) {
          seq2 <- paste0(seq2, "C")
        } else if (x > TA + TG + TC && x <= TA + TG + TC + TT) {
          seq2 <- paste0(seq2, "T")
        }
      } else {
        seq2 <- paste0(seq2, third_base)
      }
    } else {
      seq2 <- paste0(seq2, third_base)
    }
  }
  return(seq2)
}

#########################
#####   Execution   #####

# -----------------------------------------------------------
# Option Parsing
# -----------------------------------------------------------

cli::cli_rule()
cat("\n")
# Show function
cli::cli_alert(text = "Launching {.emph {.pkg cyduR}}")

option_list <- list(
  make_option(c("-i", "--input_file_name"), type = "character",
              help = "Input file name (FASTA format)"),
  make_option(c("-s", "--random_type"), type = "character", default = "n3",
              help = "Shuffling method ('n3', 'gc3', 'dn23')"),
  make_option(c("-r", "--reps"), type = "integer", default = 1000,
              help = "Replications of shuffled sequences"),
  make_option(c("-m", "--motifs"), type = "character", default = "./motif.txt",
              help = "Path to motif file"),
  make_option(c("-o", "--out_folder"), type = "character", default = "./",
              help = "Path to output folder (use full path and trailing '/')"),
  make_option(c("-d", "--delete"), type = "logical", default = FALSE,
              help = "Delete intermediate shuffled files (TRUE/FALSE)"),
  make_option(c("--seed"), type = "integer", default = 99,
              help = "Random seed for reproducibility")
)

parser <- OptionParser(option_list = option_list)
opts <- parse_args(parser)

# -----------------------------------------------------------
# Validate Inputs
# -----------------------------------------------------------

safe_run({
  cli_h1("Validations")
  
  # Check input file exists
  if (!file.exists(opts$input_file_name)) {
    stop("Input file does not exist: ", opts$input_file_name)
  }
  cli_alert_success("Input FASTA file\t[ {opts$input_file_name} ]")
  
  # Check output folder
  if (!dir.exists(opts$out_folder)) {
    dir.create(opts$out_folder, recursive = TRUE)
    cli_alert_info("Output folder created: {opts$out_folder}")
  }
  cli_alert_success("Output folder\t\t[ {opts$out_folder} ]")
  
  # Validate random type
  valid_types <- c("n3", "gc3", "dn23")
  if (!opts$random_type %in% valid_types) {
    stop("Invalid shuffling method. Use one of: ", paste(valid_types, collapse = ", "))
  }
  cli_alert_success("Shuffling algorithm\t[ {opts$random_type} ]")
  
  # Validate random seed
  if (!is.numeric(opts$seed)) {
    stop("Random seed should be numeric.")
  }
  set.seed(opts$seed)
  cli_alert_success("Random seed\t\t[ {opts$seed} ]")
  
  # Validate replicates
  if (!is.numeric(opts$reps)) {
    stop("Replicates should be numeric.")
  }
  cli_alert_success("Replicates\t\t[ {opts$reps} ]")
  
  # Test success
  cli_alert_success("All checks passed!")
}, context = "Validating inputs")

# -----------------------------------------------------------
# Main Script Logic

# Nucleotides
nts = c("A","C","G","T")

# Amino Acid Translation Table
tt <- list(
  "TTT" = "F|Phe", "TTC" = "F|Phe", "TTA" = "L|Leu", "TTG" = "L|Leu",
  "TCT" = "S|Ser", "TCC" = "S|Ser", "TCA" = "S|Ser", "TCG" = "S|Ser",
  "TAT" = "Y|Tyr", "TAC" = "Y|Tyr", "TAA" = "*|Stp", "TAG" = "*|Stp",
  "TGT" = "C|Cys", "TGC" = "C|Cys", "TGA" = "*|Stp", "TGG" = "W|Trp",
  "CTT" = "L|Leu", "CTC" = "L|Leu", "CTA" = "L|Leu", "CTG" = "L|Leu",
  "CCT" = "P|Pro", "CCC" = "P|Pro", "CCA" = "P|Pro", "CCG" = "P|Pro",
  "CAT" = "H|His", "CAC" = "H|His", "CAA" = "Q|Gln", "CAG" = "Q|Gln",
  "CGT" = "R|Arg", "CGC" = "R|Arg", "CGA" = "R|Arg", "CGG" = "R|Arg",
  "ATT" = "I|Ile", "ATC" = "I|Ile", "ATA" = "I|Ile", "ATG" = "M|Met",
  "ACT" = "T|Thr", "ACC" = "T|Thr", "ACA" = "T|Thr", "ACG" = "T|Thr",
  "AAT" = "N|Asn", "AAC" = "N|Asn", "AAA" = "K|Lys", "AAG" = "K|Lys",
  "AGT" = "S|Ser", "AGC" = "S|Ser", "AGA" = "R|Arg", "AGG" = "R|Arg",
  "GTT" = "V|Val", "GTC" = "V|Val", "GTA" = "V|Val", "GTG" = "V|Val",
  "GCT" = "A|Ala", "GCC" = "A|Ala", "GCA" = "A|Ala", "GCG" = "A|Ala",
  "GAT" = "D|Asp", "GAC" = "D|Asp", "GAA" = "E|Glu", "GAG" = "E|Glu",
  "GGT" = "G|Gly", "GGC" = "G|Gly", "GGA" = "G|Gly", "GGG" = "G|Gly"
)

safe_run({
  cli_h1("Starting Analysis")
  
  # Add the sequence ID and full string to the data 
  parsed_fasta <- read_fasta(opts$input_file_name)
  names_list <- parsed_fasta$names_list
  data <- parsed_fasta$data
  
  # Run the core workflow on all input sequences
  for(i in 1:length(names_list)) {
    # 1. - Make copy FASTA of each individual sequence in the input .fasta
    # Separate all input sequence IDs to per ID FASTA .fas files
    # Create compatible file names
    seq_name_parts <- stringr::str_split(string = names_list[i], pattern = " ", simplify = T)
    
    # Write the sequence to a .fas file
    if(ncol(seq_name_parts) == 1) {
      fas_file_name <- paste0(opts$out_folder, seq_name_parts[1,1], '.fas')
    } else if(ncol(seq_name_parts) > 1) {
      fas_file_name <- paste0(opts$out_folder, paste0(seq_name_parts[1:3], collapse = "_"), '.fas')
    }
    
    cli_alert_info("Writing [ {fas_file_name} ]")
    # Write the FASTA name and sequence to the new .fas file
    append_fasta_entry(fasta_file = fas_file_name,
                       entry = paste0(">", names_list[i], "\n", data[[eval(names_list[i])]]))
    
    # 2. - Instantiate the shuffling FASTA with the wild-type sequence
    # Write the wild-type sequence to a shuffle specific file
    shuffle_file_name <- stringr::str_replace(string = fas_file_name,
                                              pattern = ".fas",
                                              replacement = paste0("_",opts$random_type,".fasta"))
    cli_alert_info("Writing [ {shuffle_file_name} ]")
    append_fasta_entry(fasta_file = shuffle_file_name,
                       entry = paste0(">", names_list[i], "\n", data[[eval(names_list[i])]]))
    cat("\n")
    
    # 3. - Shuffle the FASTA sequences and append to shuffle FASTA
    # Run the shuffling algorithm
    # Append permuted sequences to the FASTA file
    for (j in 1:opts$reps) {
      # Start with the original sequence
      outseq <- data[[eval(names_list[i])]]
      
      # Apply permutation function based on random_type
      if (opts$random_type == "gc3") {
        outseq <- gc3(outseq)
        
      } else if(opts$random_type == "n3") {
        outseq <- n3(outseq)
        
      } else if(opts$random_type == "dn23") {
        outseq <- dn23(outseq)
      }
      
      # Output progress
      if(j %% 50 == 0) {
        cli_alert_info("[ {j} ] shuffling iterations complete ...")
      }
      # Append the permuted sequence to the file
      append_fasta_entry(fasta_file = shuffle_file_name,
                         entry = paste0(">replicate_", j, "\n", outseq, "\n"))
    }
    
    # 4. - Run the statistics reporter
    # Running the stats testing
    results_file_name <- stringr::str_replace(string = shuffle_file_name,
                                              pattern = ".fasta",
                                              replacement = "_results.txt")
    if(!is.null(opts$motifs)) {
      # Build the command to run shmsim with motifs provided
      cmd <- paste("~/CDUR/shmsim", 
                   shuffle_file_name,
                   opts$motifs, 
                   ">", 
                   results_file_name)
      
      # Run the command
      cat("\n")
      cli_alert_info("Running command:\n\n\t{cmd}")
      #system(cmd, intern = FALSE)
      
      # Check if delete flag is set to TRUE
      # run the clean-up
      if (opts$delete) {
        # Deleting the .fas file
        if (file.exists(fas_file_name)) {
          file.remove(fas_file_name)
        }
        # deleting the .fasta file
        if (file.exists(shuffle_file_name)) {
          file.remove(shuffle_file_name)
        }
      }
    } else {
      # If motifs are not provided, run shmsim without motifs
      cmd <- paste("~/CDUR/shmsim", 
                   shuffle_file_name,
                   ">",
                   results_file_name)
      
      # Run the command
      cat("\n")
      cli_alert_info("Running command:\n\n\t{cmd}")
      #system(cmd, intern = FALSE)
    }
  }
  
  cli_alert_success("Analysis completed!")
}, context = "Running main script")

quit(save = "no", status = 0)
