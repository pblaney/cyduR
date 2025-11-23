#!/usr/bin/env bash
# Install necessary R packages

echo "Installing necessary R packages ..."
echo 

cran_deps=(optparse cli clisymbols stringr)
for i in "${cran_deps[@]}"
    do
        echo 
        echo "INSTALL - [ ${i} ]"
        cran_cmd="R -q -e 'remotes::install_cran(\"${i}\")'"
        (eval $cran_cmd)
        sleep 2
    done

echo 
echo "Installation of R dependencies complete ..."
