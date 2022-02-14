#!/bin/bash


for TRP in {8..12}; do
    python Demo_STP.py $TRP 0
    python Demo_STP.py $TRP 1
done
