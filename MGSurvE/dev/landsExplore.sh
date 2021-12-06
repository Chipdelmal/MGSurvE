#!/bin/bash

N=2
PTH_O='/RAID5/marshallShare/MGS_Demos'
LND='UNIF'

LG='\033[1;34m'
NC='\033[0m'

(
for id in {50..54}; do
    ((i=i%N)); ((i++==0)) && wait
    echo -e "${LG}* Launched landscape: ${id} ${NC}"
    python Optim_dev.py  $PTH_O $LND $id &
done
)
