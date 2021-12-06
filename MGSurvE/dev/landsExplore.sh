#!/bin/bash

PTH_O='/RAID5/marshallShare/MGS_Demos'
LND='UNIF'

LG='\033[1;34m'
NC='\033[0m' # No Color

for id in {1..4}
do
    echo -e "${LG}* Launched landscape: ${id} ${NC}"
    python Optim_dev.py  $PTH_O $LND $id
    echo -e "${LG}* Finished landscape: ${id} ${NC}"
done
