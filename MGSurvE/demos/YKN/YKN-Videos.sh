#!/bin/bash

WHT='\033[0;37m'  
RED='\033[0;31m'
NCL='\033[0m'


for i in {1..10}; do
    printf "${RED}* ID $i ...${NCL}\n" 
    python YKN-Video.py YKND man "$i"
done

for i in {1..10}; do
    printf "${RED}* ID $i ...${NCL}\n" 
    python YKN-Video.py YKND max "$i"
done

for i in {1..10}; do
    printf "${RED}* ID $i ...${NCL}\n" 
    python YKN-Video.py YKNC man "$i"
done

for i in {1..10}; do
    printf "${RED}* ID $i ...${NCL}\n" 
    python YKN-Video.py YKNC max "$i"
done