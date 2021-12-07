#!/bin/bash

N=2
PTH_O='/RAID5/marshallShare/MGS_Demos'
LND='DNUT'

LG='\033[1;34m'
NC='\033[0m'

# (
for id in {3..6}; do
    # ((i=i%N)); ((i++==0)) && wait
    ID="S${id}"
    echo -e "${LG}* Launched landscape: ${ID} ${NC}"
    python Optim_dev.py  $PTH_O $LND $ID 
done
# )
