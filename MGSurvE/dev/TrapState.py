#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math
import numpy as np
import pandas as pd
from os import path
from sys import argv
from copy import deepcopy
import matplotlib.pyplot as plt
from deap import base, creator, algorithms, tools
from compress_pickle import dump, load
import MGSurvE as srv
import warnings
warnings.filterwarnings('ignore', 'The iteration is not making good progress')


(OUT_PTH, LND_TYPE, ID, PTS_TYPE) = ('./Lands', 'AGG', '001', 3)

ptsNum = 20
bbox = ((-225, 225), (-175, 175))
xy = srv.ptsRandUniform(ptsNum, bbox).T
pType = np.random.choice(PTS_TYPE, xy.shape[1])
points = pd.DataFrame({'x': xy[0], 'y': xy[1], 't': pType})
###############################################################################
# Defining Traps
###############################################################################
traps = pd.DataFrame({
    'x': [-1, 0, 2, 0],
    'y': [0, 1, 0, -1],
    't': [0, 1, 0, 1],
    'f': [0, 0, 0, 0]
})
tker = {
    0: {'kernel': srv.exponentialDecay, 'params': {'A': 0.50, 'b': .050}},
    1: {'kernel': srv.exponentialDecay, 'params': {'A': 0.35, 'b': .025}}
}
tMask = np.asarray([
    [.9, .1, 0],
    [0, .8, .2]
])
###############################################################################
# Landscape
###############################################################################
lnd = srv.Landscape(
    points, 
    kernelParams={'params': [.075, 1.0e-10, math.inf], 'zeroInflation': .75},
    traps=traps, trapsKernels=tker
)
###############################################################################
# Masked trapping
###############################################################################
pointTypes = lnd.pointTypes
trapsTypes = lnd.trapsTypes
trapsKernels = lnd.trapsKernels
trapsDistances = lnd.trapsDistances


# Base unbiased probs ---------------------------------------------------------
trapProbs = np.asarray([
    [
        trapsKernels[ttype]['kernel'](i, **trapsKernels[ttype]['params'])
        for (i, ttype) in zip(dist, trapsTypes)
    ] for dist in trapsDistances
])
# Point-type to trap probs ----------------------------------------------------
trapMask = np.asarray([
    [
        tMask[ttype][pointTypes[ix]]
        for ttype in trapsTypes
    ] for ix in range(len(trapsDistances))
])
# Alpha -----------------------------------------------------------------------
trapMatrix = trapProbs*trapMask