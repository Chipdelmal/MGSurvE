#!/bin/bash

WHT='\033[0;37m'  
RED='\033[0;31m'
NCL='\033[0m'


for trp in 10 15 20; do
    printf "* TRPS $trp ...\n" 
    for i in {1..10}; do
        printf "${RED}\t* ID $i ...${NCL}\n" 
        python STP-Video.py STPD man "$trp" "$i"
    done
done

