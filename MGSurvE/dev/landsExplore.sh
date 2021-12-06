#!/bin/bash

PTH_O='/RAID5/marshallShare/MGS_Demos'
LND='UNIF'

for id in {1..5}
do
    python Optim_dev.py  $PTH_O $LND $id
    echo "Finished landscape: ${id}"
done