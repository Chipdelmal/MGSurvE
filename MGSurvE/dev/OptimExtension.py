#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math
import random
import numpy as np
import pandas as pd
from os import path
from sys import argv
from copy import deepcopy
from deap import tools
import matplotlib.pyplot as plt
import MGSurvE as srv
import numpy as np


(OUT_PTH, LND_TYPE, ID) = ('./Lands', 'DNUT', 'D01')
###############################################################################
# Defining Landscape and Traps
###############################################################################
if LND_TYPE == 'UNIF':
    ptsNum = 100
    bbox = ((-225, 225), (-175, 175))
    xy = srv.ptsRandUniform(ptsNum, bbox).T
elif LND_TYPE == 'GRID':
    ptsNum = 20
    bbox = ((-225, 225), (-225, 225))
    xy = srv.ptsRegularGrid(ptsNum, bbox).T
elif LND_TYPE == 'DNUT':
    ptsNum = 50
    radii = (100, 150)
    xy = srv.ptsDonut(ptsNum, radii).T
points = pd.DataFrame({'x': xy[0], 'y': xy[1], 't': [0]*xy.shape[1]})
# Traps info ------------------------------------------------------------------
traps = pd.DataFrame({
    'x': [0, 0, 0, 0, 0],
    'y': [0, 0, 0, 0, 0],
    't': [1, 3, 2, 1, 0],
    'f': [0, 0, 1, 0, 0],
    'o': [1, 1, 0, 0, 0]
})
tKernels = {
    0: {'kernel': srv.exponentialDecay, 'params': {'A': .3, 'b': .05}},
    1: {'kernel': srv.exponentialDecay, 'params': {'A': .35, 'b': .04}},
    2: {'kernel': srv.exponentialDecay, 'params': {'A': .25,  'b': .025}} ,
    3: {'kernel': srv.sigmoidDecay,     'params': {'A': .2, 'rate': .5, 'x0': 1}}
}
###############################################################################
# Defining Landscape and Traps
###############################################################################
lnd = srv.Landscape(
    points, 
    kernelParams={'params': [.075, 1.0e-10, math.inf], 'zeroInflation': .75},
    traps=traps, trapsKernels=tKernels
)
bbox = lnd.getBoundingBox()
trpMsk = srv.genFixedTrapsMask(lnd.trapsFixed)
trpTsk = lnd.trapsTOptim
###############################################################################
# Optimization Extension
###############################################################################
lndGA = deepcopy(lnd)
chromosome = srv.initChromosomeMixed(
    trapsCoords=lndGA.trapsCoords, fixedTrapsMask=trpMsk, coordsRange=bbox,
    trapsPool=[0, 1, 2, 3, 0, 1, 2]   
)
# Develop mutation operator ---------------------------------------------------
(coordSect, typesSect) = (chromosome[:len(trpMsk)], chromosome[len(trpMsk):])
coordSect = srv.mutateChromosome(coordSect, trpMsk)

typesSect
mutShuffleIndexes(typesSect, 1, trpTsk)

(coordSect[0]+typesSect[0], )





len(typesSect)