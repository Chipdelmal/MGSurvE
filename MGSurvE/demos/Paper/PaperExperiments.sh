#!/bin/bash

declare -a lnds=("Grid" "Uniform" "Ring" "Poisson" )
###############################################################################
# Setting landscapes up
###############################################################################
echo "  * [1/2] Generating landscapes..."
for lnd in ${lnds[@]}; do
    python Landscape.py $lnd
done
###############################################################################
# Optimizing traps
###############################################################################
echo "  * [2/2] Optimizing landscape..."
for lnd in ${lnds[@]}; do
    python Optimization.py "${lnd}_LND_HOM"
    python Optimization.py "${lnd}_LND_HET"
done

