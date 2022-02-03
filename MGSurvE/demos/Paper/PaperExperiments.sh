#!/bin/bash

###############################################################################
# Setting landscapes up
###############################################################################
echo "  * [1/3] Generating landscapes..."
python Grid_Landscape.py
python Ring_Landscape.py
python Uniform_Landscape.py
python Poisson_Landscape.py
###############################################################################
# Optimizing traps
###############################################################################
echo "  * [2/3] Optimizing landscape (1/2)..."
python Optimization.py 'Grid_LND_HOM'
python Optimization.py 'Ring_LND_HOM'
python Optimization.py 'Uniform_LND_HOM'
python Optimization.py 'Poisson_LND_HOM'
echo "  * [2/3] Optimizing landscape (2/2)..."
python Optimization.py 'Grid_LND_HET'
python Optimization.py 'Ring_LND_HET'
python Optimization.py 'Uniform_LND_HET'
python Optimization.py 'Poisson_LND_HET'