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
    eval(expr),
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

#########################
#####   Execution   #####

# -----------------------------------------------------------
# Option Parsing
# -----------------------------------------------------------

option_list <- list(
  make_option(c("-i", "--input_file_name"), type = "character", required = TRUE,
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
  cli_alert_success("Input FASTA file [ {opts$input_file_name} ]")
  
  # Check output folder
  if (!dir.exists(opts$out_folder)) {
    dir.create(opts$out_folder, recursive = TRUE)
    cli_alert_info("Output folder created: {opts$out_folder}")
  }
  cli_alert_success("Output folder [ {opts$out_folder} ]")
  
  # Validate random type
  valid_types <- c("n3", "gc3", "dn23")
  if (!opts$random_type %in% valid_types) {
    stop("Invalid shuffling method. Use one of: ", paste(valid_types, collapse = ", "))
  }
  cli_alert_success("Shuffling algorithm [ {opts$random_type} ]")
  
  # Validate random seed
  if (!is.numeric(opts$seed)) {
    stop("Random seed should be numeric.")
  }
  set.seed(opts$seed)
  cli_alert_success("Random seed [ {opts$seed} ]")
  
  # Test success
  cli_alert_success("All checks passed!")
}, context = "Validating inputs")

# -----------------------------------------------------------
# Main Script Logic

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
    if(ncol(seq_name_parts) < 4) {
      fas_file_name <- paste0(opts$out_folder, paste0(seq_name_parts, collapse = "_"), '.fas')
    } else {
      fas_file_name <- paste0(opts$out_folder, paste0(seq_name_parts[1:4], collapse = "_"), '.fas')
    }
    cli_alert_info("Writing [ ", fas_file_name, " ]")
    # Write the FASTA name and sequence to the new .fas file
    append_fasta_entry(fasta_file = fas_file_name,
                       entry = paste0(">", names_list[i], "\n", data[[eval(names_list[i])]]))
    
    # 2. - Instantiate the shuffling FASTA with the wild-type sequence
    # Write the wild-type sequence to a shuffle specific file
    shuffle_file_name <- stringr::str_replace(string = fas_file_name,
                                              pattern = ".fas",
                                              replacement = paste0("_",opts$random_type,".fasta"))
    cli_alert_info("Writing [ ", shuffle_file_name, " ]")
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
        outseq <- paste0("dn23_TESTING", j)
        #outseq <- dn23(outseq)
      }
      
      # Output progress
      if(j %% 50 == 0) {
        cli_alert_info("[ ", j, " ] shuffling iterations complete ...")
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
      cli_alert_info("Running command:\n  ",cmd)
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
      cli_alert_info("Running command:\n  ",cmd)
      #system(cmd, intern = FALSE)
    }
  }
  
  cli_alert_success("Analysis completed!")
}, context = "Running main script")

quit(save = "no", status = 0)
