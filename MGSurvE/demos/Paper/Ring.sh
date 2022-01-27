#!/bin/bash

###############################################################################
# Setting landscapes up
###############################################################################
echo "  * [1/x] Generating landscapes..."
python Ring_Landscape.py
###############################################################################
# Optimizing traps
###############################################################################
echo "  * [2/x] Optimizing landscape (1/2)..."
python Ring_Optimization.py 'Ring_LND_HOM'
echo "  * [2/x] Optimizing landscape (2/2)..."
python Ring_Optimization.py 'Ring_LND_HET'