#!/bin/bash

# scp -r lab:/RAID5/marshallShare/MGS_Benchmarks/Grid/* /home/chipdelmal/Documents/WorkSims/MGSurvE_Benchmarks/Grid

# PTH_O='/home/chipdelmal/Documents/WorkSims/MGSurvE_Benchmarks'
PTH_O=/RAID5/marshallShare/MGS_Benchmarks/Grid
GENS=1000


LG='\033[1;34m'
NC='\033[0m'


for pts in 50 100 200 300 400; do
    for trps in 1 2 5 10 15 20; do
        for rep in {1..10}; do
            echo -e "${LG}* LND: ${trps} ${pts} ${rep} ${NC}"
            python GridOptimization.py $trps $pts $rep $GENS $PTH_O
        done
    done
done
