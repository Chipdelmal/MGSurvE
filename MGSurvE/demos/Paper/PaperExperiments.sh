#!/bin/bash

OPT=$1
declare -a lnds=("Grid" "Uniform" "Ring" "Poisson" )
###############################################################################
# Setting landscapes up
###############################################################################
echo "* [1/3] Generating landscapes"
# for lnd in ${lnds[@]}; do
#     printf "\r\tGenerating $lnd..."
#     python Landscape.py $lnd "ZI"
#     python Landscape.py $lnd "ZN"
#     printf "\r\033[K"
# done
###############################################################################
# Optimizing Continuous Traps
###############################################################################
echo "* [2/3] Optimizing landscapes"
if [ $OPT == "Simple" ];
then
    for lnd in ${lnds[@]}; do
        printf "\r\tOptimizing (Simple DO) $lnd Zero-Inflated..."
        python Optimization-Simple.py "${lnd}_LND_HOM" "ZI"
        python Optimization-Simple.py "${lnd}_LND_HET" "ZI"
        printf "\r\tOptimizing (Simple DO) $lnd Non Zero-Inflated..."
        python Optimization-Simple.py "${lnd}_LND_HOM" "ZN"
        python Optimization-Simple.py "${lnd}_LND_HET" "ZN"
        printf "\r\033[K"
    done
else
    for lnd in ${lnds[@]}; do
        printf "\r\tOptimizing (Simple CO) $lnd Zero-Inflated..."
        python Optimization.py "${lnd}_LND_HOM" "ZI"
        python Optimization.py "${lnd}_LND_HET" "ZI"
        printf "\r\tOptimizing (Simple CO) $lnd Non Zero-Inflated..."
        python Optimization.py "${lnd}_LND_HOM" "ZN"
        python Optimization.py "${lnd}_LND_HET" "ZN"
        printf "\r\033[K"
    done
fi
###############################################################################
# Optimizing Discrete Traps
###############################################################################
echo "* [3/3] Optimizing landscapes"
if [ $OPT == "Simple" ];
then
    for lnd in ${lnds[@]}; do
        printf "\r\tOptimizing (Complex DO) $lnd Zero-Inflated..."
        python OptimizationDO-Simple.py "${lnd}_LND_HOM" "ZI"
        python OptimizationDO-Simple.py "${lnd}_LND_HET" "ZI"
        printf "\r\tOptimizing (Complex DO) $lnd Non Zero-Inflated..."
        python OptimizationDO-Simple.py "${lnd}_LND_HOM" "ZI"
        python OptimizationDO-Simple.py "${lnd}_LND_HET" "ZI"
        printf "\r\033[K"
    done
else
    for lnd in ${lnds[@]}; do
        printf "\r\tOptimizing (Complex CO) $lnd Zero-Inflated..."
        python OptimizationDO.py "${lnd}_LND_HOM" "ZI"
        python OptimizationDO.py "${lnd}_LND_HET" "ZI"
        printf "\r\tOptimizing (Complex CO) $lnd Non Zero-Inflated..."
        python OptimizationDO.py "${lnd}_LND_HOM" "ZN"
        python OptimizationDO.py "${lnd}_LND_HET" "ZN"
        printf "\r\033[K"
    done
fi
###############################################################################
# Goodbye
###############################################################################
echo "* Done!"