#!/bin/bash

WHT='\033[0;37m'  
RED='\033[0;31m'
NCL='\033[0m'

for i in {1..30}; do
    printf "${WHT}* Iteration $i...${NCL}\n" 
    for m in "sum" "man"; do  # "max"; do 
        printf "${RED}\t* Optimizing $m... ${NCL}\n"
        python YKN-Discrete.py "YKN" "$m" "$i"
    done
done
