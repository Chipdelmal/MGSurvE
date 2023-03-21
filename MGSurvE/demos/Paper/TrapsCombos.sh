#!/bin/bash

declare -a combo=("sum-men" "men-max" "max-men" "max-max" "men-men" "men-sum" "sum-max" "max-sum")
for cmb in ${combo[@]}; do
    python OptimizationDO-Validation.py "Grid_LND_HOM" "ZI" "$cmb"
done