#!/usr/bin/env bash
# Install necessary R packages

echo "Installing necessary R packages ..."
echo 

cran_deps=(optparse cli clisymbols stringr dplyr)
for i in "${cran_deps[@]}"
    do
        echo 
        echo "INSTALL - [ ${i} ]"
        cran_cmd="R -q -e 'remotes::install_cran(\"${i}\")'"
        (eval $cran_cmd)
        sleep 2
    done

bioc_deps=(GenomeInfoDb GenomicRanges rtracklayer BSgenome.Hsapiens.UCSC.hg38 nullranges Biostrings)
for i in "${bioc_deps[@]}"
    do
        echo 
        echo "INSTALL - [ ${i} ]"
        bioc_cmd="R -q -e 'remotes::install_bioc(\"${i}\")'"
        (eval $bioc_cmd)
        sleep 2
    done

echo 
echo "Installation of R dependencies complete ..."
