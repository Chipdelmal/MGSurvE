#!/bin/bash

declare -a lnds=("Grid" "Uniform" "Ring" "Poisson" )
###############################################################################
# Setting landscapes up
###############################################################################
echo "* [1/2] Generating landscapes"
for lnd in ${lnds[@]}; do
    printf "\r\tGenerating $lnd..."
    python Landscape.py $lnd
    printf "\r\033[K"
done
###############################################################################
# Optimizing traps
###############################################################################
echo "* [2/2] Optimizing landscapes"
for lnd in ${lnds[@]}; do
    printf "\r\tOptimizing $lnd..."
    python Optimization.py "${lnd}_LND_HOM"
    python Optimization.py "${lnd}_LND_HET"
    printf "\r\033[K"
done
