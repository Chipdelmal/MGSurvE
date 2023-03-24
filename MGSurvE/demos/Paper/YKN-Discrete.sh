#!/bin/bash

WHT='\033[0;37m'  
RED='\033[0;31m'
NCL='\033[0m'

for m in "sum" "man" "max"; do
    printf "${WHT}* Optimizing $m...${NCL}\n"
    for i in {1..30}; do
        printf "${RED}\t* Running iteration $i... ${NCL}\n"
        python YKN-Discrete.py "YKN" "$m" "$i"
    done
done
