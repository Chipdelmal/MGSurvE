#!/bin/bash

OPT=$1
declare -a lnds=("Grid" "Uniform" "Ring" "Poisson" )
###############################################################################
# Setting landscapes up
###############################################################################
echo "* [1/3] Generating landscapes"
for lnd in ${lnds[@]}; do
    printf "\r\tGenerating $lnd..."
    python Landscape.py $lnd
    printf "\r\033[K"
done
###############################################################################
# Optimizing Continuous Traps
###############################################################################
echo "* [2/3] Optimizing landscapes"
# if [ $OPT == "Simple" ];
# then
    for lnd in ${lnds[@]}; do
        printf "\r\tOptimizing (Simple DO) $lnd..."
        python OptimizationDO-Simple.py "${lnd}_LND_HOM"
        python OptimizationDO-Simple.py "${lnd}_LND_HET"
        printf "\r\033[K"
    done
# else
    for lnd in ${lnds[@]}; do
        printf "\r\tOptimizing (Simple CO) $lnd..."
        python Optimization-Simple.py "${lnd}_LND_HOM"
        python Optimization-Simple.py "${lnd}_LND_HET"
        printf "\r\033[K"
    done
# fi
###############################################################################
# Optimizing Discrete Traps
###############################################################################
echo "* [3/3] Optimizing landscapes"
# if [ $OPT == "Simple" ];
# then
    for lnd in ${lnds[@]}; do
        printf "\r\tOptimizing (Complex DO) $lnd..."
        python OptimizationDO.py "${lnd}_LND_HOM"
        python OptimizationDO.py "${lnd}_LND_HET"
        printf "\r\033[K"
    done
# else
    for lnd in ${lnds[@]}; do
        printf "\r\tOptimizing (Complex CO) $lnd..."
        python Optimization.py "${lnd}_LND_HOM"
        python Optimization.py "${lnd}_LND_HET"
        printf "\r\033[K"
    done
# fi
###############################################################################
# Goodbye
###############################################################################
echo "* Done!"