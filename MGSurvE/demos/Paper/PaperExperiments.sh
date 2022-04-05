#!/bin/bash

OPT=$1
declare -a lnds=("Grid" "Uniform" "Ring" "Poisson" )
###############################################################################
# Setting landscapes up
###############################################################################
# echo "* [1/3] Generating landscapes"
# for lnd in ${lnds[@]}; do
#     printf "\r\tGenerating $lnd..."
#     python Landscape.py $lnd
#     printf "\r\033[K"
# done
###############################################################################
# Optimizing traps
###############################################################################
echo "* [2/3] Optimizing landscapes"
if [ $OPT == "Simple" ];
then
    for lnd in ${lnds[@]}; do
        printf "\r\tOptimizing (Simple) $lnd..."
        python Optimization.py "${lnd}_LND_HOM"
        python Optimization.py "${lnd}_LND_HET"
        printf "\r\033[K"
    done
else
    for lnd in ${lnds[@]}; do
        printf "\r\tOptimizing (Complex) $lnd..."
        python Optimization.py "${lnd}_LND_HOM"
        python Optimization.py "${lnd}_LND_HET"
        printf "\r\033[K"
    done
fi
###############################################################################
# Concatenating results
###############################################################################
echo "* [3/3] Concatenating images"
for lnd in ${lnds[@]}; do
    printf "\r\tConcatenating $lnd..."
    convert "./sims_out/${lnd}_LND_HOM_TRP.png" "./sims_out/${lnd}_LND_HET_TRP.png" +append "./sims_out/${lnd}_LND_TRP.png"
    printf "\r\033[K"
done
