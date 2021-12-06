#!/bin/bash

PTH_O='/RAID5/marshallShare/MGS_Demos'
LND='UNIF'


N=2
(
for id in {1..4}; do
    ((i=i%N)); ((i++==0)) && wait
    python Optim_dev.py  $PTH_O $LND $id &
    echo "Finished landscape: ${id}"
done
)