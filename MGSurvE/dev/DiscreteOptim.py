#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math
from operator import index
import numpy as np
import pandas as pd
from os import path
from sys import argv
import random
from random import choice
import matplotlib.pyplot as plt
import MGSurvE as srv
import numpy as np
import warnings
warnings.filterwarnings('ignore', 'The iteration is not making good progress')

if srv.isNotebook():
    (OUT_PTH, LND_TYPE, ID) = ('./Lands', 'GRID', 'G01')
else:
    (OUT_PTH, LND_TYPE, ID) = (argv[1], argv[2], argv[3].zfill(3))
###############################################################################
# Defining Landscape and Traps
###############################################################################
if LND_TYPE == 'UNIF':
    ptsNum = 400
    bbox = ((-225, 225), (-175, 175))
    xy = srv.ptsRandUniform(ptsNum, bbox).T
elif LND_TYPE == 'GRID':
    ptsNum = 10
    bbox = ((-225, 225), (-225, 225))
    xy = srv.ptsRegularGrid(ptsNum, bbox).T
elif LND_TYPE == 'DNUT':
    ptsNum = 150
    radii = (100, 150)
    xy = srv.ptsDonut(ptsNum, radii).T
points = pd.DataFrame({
    'x': xy[0], 'y': xy[1], 
    't': [0]*xy.shape[1], 'id': range(0, 100)
})
# Traps info ------------------------------------------------------------------
traps = pd.DataFrame({
    'sid': [0, 10, 55, 25],
    'x': [xy[0, 0], 0, 0, 0],
    'y': [xy[1, 0], 0, 0, 0],
    't': [0, 0, 0 ,0],
    'f': [0, 0, 0 ,1]
})
tKernels = {0: {'kernel': srv.exponentialDecay, 'params': {'A': .2, 'b': 0.1}}}
###############################################################################
# Defining Landscape and Traps
###############################################################################
lnd = srv.Landscape(
    points, kernelParams={'params': [.075, 1.0e-10, math.inf], 'zeroInflation': .75},
    traps=traps, trapsKernels=tKernels, pointsTrapBanned={5}
)
srv.dumpLandscape(lnd, OUT_PTH, '{}_{}_CLN'.format(LND_TYPE, ID))
bbox = lnd.getBoundingBox()
trpMsk = srv.genFixedTrapsMask(lnd.trapsFixed)
###############################################################################
# Implementing extension
###############################################################################
lnd.calcFundamentalMatrix()
lnd.getDaysTillTrapped()
# lnd.pointCoords
lnd.fundamentalMatrix
# Init chromosome -------------------------------------------------------------
trpsIDPos  = [0, 10, 55, 25]
fixedTraps = [0, 0, 0, 1]
trapsNum = lnd.trapsNumber
ptsNum = lnd.pointNumber
ptsIds = tuple((range(ptsNum)))

def calcDiscreteFitness(
        chromosome, landscape,
        optimFunction=srv.getDaysTillTrapped,
        optimFunctionArgs={'outer': np.mean, 'inner': np.max},
    ):
    ptsIds = landscape.pointID
    siteIndex = [ptsIds.index(i) for i in chromosome]
    trapXY = np.asarray([landscape.pointCoords[i] for i in siteIndex])
    fit = srv.calcFitness(
        trapXY, landscape=landscape,
        optimFunction=optimFunction,
        optimFunctionArgs=optimFunctionArgs
    )
    return fit

def calcDiscreteFitnessPseudoInverse(
        chromosome, landscape,
        optimFunction=srv.getDaysTillTrappedPseudoInverse,
        optimFunctionArgs={'outer': np.mean, 'inner': np.max},  
        rcond=1e-30
    ):
    ptsIds = landscape.pointID
    siteIndex = [ptsIds.index(i) for i in chromosome]
    trapXY = np.asarray([landscape.pointCoords[i] for i in siteIndex])
    fit = srv.calcFitnessPseudoInverse(
        trapXY, landscape=landscape,
        optimFunction=optimFunction,
        optimFunctionArgs=optimFunctionArgs,
        rcond=rcond
    )
    return fit


chromB = srv.initDiscreteChromosome(lnd.pointID, lnd.trapsFixed, lnd.pointsTrapBanned)
chromA = srv.mutateDiscreteChromosome(
    chromB.copy(), lnd.pointID, lnd.trapsFixed, indpb=1
)[0]
print(chromA, chromB)
print(cxUniform(chromA, chromB,  lnd.trapsFixed, indpb=.5))
calcDiscreteFitness(chromA, lnd)
calcDiscreteFitnessPseudoInverse(chromB, lnd)

chromosome = chromA
landscape = lnd

siteIndex = [ptsIds.index(i) for i in chromosome]
trapXY = np.asarray([landscape.pointCoords[i] for i in siteIndex])
np.reshape(trapXY, (-1, 2))


