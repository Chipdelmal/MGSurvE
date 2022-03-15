#!/bin/bash


for TRP in {10..20..2}; do
    python Demo_STP.py $TRP 0
    python Demo_STP.py $TRP 1
done
