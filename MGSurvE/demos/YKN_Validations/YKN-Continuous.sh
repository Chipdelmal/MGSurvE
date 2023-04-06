#!/bin/bash

WHT='\033[0;37m'  
RED='\033[0;31m'
NCL='\033[0m'

for i in {1..10}; do
    printf "${WHT}* Iteration $i...${NCL}\n" 
    for m in "man"; do # "max" "sum"; do 
        printf "${RED}\t* Optimizing $m... ${NCL}\n"
        python YKN-Continuous.py "YKN" "$m" "$i"
    done
done
