#!/bin/bash

WHT='\033[0;37m'  
RED='\033[0;31m'
NCL='\033[0m'

for trps in 5 10 15 20; do
    for id in {1..10}; do
        printf "${WHT}* Iteration $i...${NCL}\n" 
        for m in "man"; do 
            printf "${RED}\t* Optimizing $m... ${NCL}\n"
            python STP-Discrete.py "$m" "$trps" "$id"
        done
    done
done
