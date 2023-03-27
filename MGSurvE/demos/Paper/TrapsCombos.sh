#!/bin/bash

WHT='\033[0;37m'  
RED='\033[0;31m'
NCL='\033[0m'

declare -a combo=('sum-men' 'sum-max' 'sum-sum' 'sum-min' 'men-men' 'men-max' 'men-sum' 'men-min')
for i in {1..30}; do
    printf "${WHT}* Running iteration $i... ${NCL}\n"
    for cmb in ${combo[@]}; do
        printf "${RED}\t* Optimizing $cmb...${NCL}\n"
        python OptimizationDO-Validation.py "Grid_LND_HOM" "ZI" "$cmb" "$i"
    done
done