#!/bin/bash

declare -a combo=('sum-men' 'sum-max' 'sum-sum' 'sum-min' 'men-men' 'men-max' 'men-sum' 'men-min')
for cmb in ${combo[@]}; do
    python OptimizationDO-Validation.py "Grid_LND_HOM" "ZI" "$cmb"
done
