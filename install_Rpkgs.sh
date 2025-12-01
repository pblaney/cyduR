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

echo "INSTALL - [ BSgenome.Hsapiens.UCSC.hg38 ]"
bsgenome_cmd="R -q -e 'BiocManager::install(\"BSgenome.Hsapiens.UCSC.hg38\")'"
(eval $bsgenome_cmd)
sleep 2

echo 
echo "INSTALL - [ nullranges ]"
nullranges_cmd="R -q -e 'remotes::install_bioc(\"nullranges\")'"
(eval $nullranges_cmd)
sleep 2

echo 
echo "Installation of R dependencies complete ..."
