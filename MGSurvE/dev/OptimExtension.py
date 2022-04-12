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
import numpy.random as rand


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
    'x': [0, 0, 0, 0, 0, 0],
    'y': [0, 0, 0, 0, 0, 0],
    't': [0, 1, 2, 3, 2, 1],
    'f': [1, 1, 1, 1, 1, 1],
    'o': [0, 0, 0, 0, 0, 0]
})
tKer = {
    0: {'kernel': srv.exponentialDecay, 'params': {'A': .3, 'b': .05}},
    1: {'kernel': srv.exponentialDecay, 'params': {'A': .35, 'b': .04}},
    2: {'kernel': srv.exponentialDecay, 'params': {'A': .25,  'b': .025}} ,
    3: {'kernel': srv.sigmoidDecay,     'params': {'A': .2, 'rate': .5, 'x0': 1}}
}
###############################################################################
# Defining Landscape and Traps
###############################################################################
lnd = srv.Landscape(points, traps=traps, trapsKernels=tKer)
lndBase = deepcopy(lnd)
bbox = lnd.getBoundingBox()
###############################################################################
# Optimization Extension
###############################################################################
trapsPool = list(traps['t'])+list(range(100, 110))
trpsNum = traps.shape[0]
baseChrom = srv.initChromosomeMixed(
    trapsCoords=lndBase.trapsCoords, 
    fixedTrapsMask=srv.genFixedTrapsMask(lndBase.trapsFixed), 
    typeOptimMask=lnd.trapsTOptim,
    coordsRange=bbox, indpb=1,
    trapsPool=trapsPool
)
# Test init routine -----------------------------------------------------------
lndTest = deepcopy(lndBase)
traps = pd.DataFrame({
    'x': [0, 0, 0, 0, 0, 0],
    'y': [0, 0, 0, 0, 0, 0],
    't': [0, 1, 2, 3, 2, 1],
    'f': [0, 1, 0, 1, 0, 1],
    'o': [1, 0, 1, 0, 1, 0]
})
lndTest.updateTraps(traps, tKer)
testChrom = srv.initChromosomeMixed(
    trapsCoords=lndTest.trapsCoords, 
    fixedTrapsMask=srv.genFixedTrapsMask(lndTest.trapsFixed), 
    typeOptimMask=list(traps['o']),
    coordsRange=bbox, indpb=1,
    trapsPool=trapsPool
)
###############################################################################
# Dev Extensions
###############################################################################
trapsCoords = lndTest.trapsCoords
fixedTrapsMask = srv.genFixedTrapsMask(lndTest.trapsFixed)
typeOptimMask = lndTest.trapsTOptim
coordsRange = bbox
indpb = .75
 
chromosome = deepcopy(baseChrom)
srv.mutateChromosomeMixed(
    chromosome, fixedTrapsMask, typeOptimMask,
    mutCoordArgs={
        'randFun': rand.normal, 
        'randArgs': {'loc': 0, 'scale': 10}, 
        'indpb': 1
    },
    mutTypeArgs={
        'indpb': 1
    }
)
