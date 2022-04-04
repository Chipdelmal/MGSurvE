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
import MGSurvE as srv
import warnings
warnings.filterwarnings('ignore', 'The iteration is not making good progress')


(OUT_PTH, LND_TYPE, ID, TRPS_NUM) = ('./scratch/', 'UNIF', 'SX01', 4)
###############################################################################
# Generating Pointsets
###############################################################################
ptsNum = 200
bbox = ((-100, 100), (-80, 80))
xy = srv.ptsRandUniform(ptsNum, bbox).T
points = pd.DataFrame({'x': xy[0], 'y': xy[1], 't': [0]*xy.shape[1]})
###############################################################################
# Defining Movement
###############################################################################
movementKernel = {
    'Male': {
        'kernelFunction': srv.zeroInflatedExponentialKernel,
        'kernelParams': {
            'params': [.050, 1.0e-10, math.inf], 'zeroInflation': .5
        }
    },
    'Female': {
        'kernelFunction': srv.zeroInflatedExponentialKernel,
        'kernelParams': {
            'params': [.025, 1.0e-10, math.inf], 'zeroInflation': .7
        }
    }
}
###############################################################################
# Defining Traps
###############################################################################
nullTraps = [0] * TRPS_NUM
traps = pd.DataFrame({
    'x': nullTraps, 'y': nullTraps, 'f': nullTraps, 't': nullTraps,
})
tKernels = {
    'Male': {
        0: {'kernel': srv.exponentialDecay, 'params': {'A': .3, 'b': .1}}
    },
    'Female': {
        0: {'kernel': srv.exponentialDecay, 'params': {'A': .75, 'b': .1}}
    }
}
###############################################################################
# Setting Landscape Up
###############################################################################
lndM = srv.Landscape(
    points, traps=traps,
    kernelFunction=movementKernel['Male']['kernelFunction'],
    kernelParams=movementKernel['Male']['kernelParams'],
    trapsKernels=tKernels['Male'], trapsRadii=[.1, ]
)
lndF = srv.Landscape(
    points, traps=traps,
    kernelFunction=movementKernel['Female']['kernelFunction'],
    kernelParams=movementKernel['Female']['kernelParams'],
    trapsKernels=tKernels['Female'], trapsRadii=[.1, ]
)
###############################################################################
# Plot Landscape
###############################################################################
(fig, ax) = plt.subplots(1, 1, figsize=(15, 15), sharey=False)
lndM.plotSites(fig, ax, size=100)
lndM.plotMigrationNetwork(fig, ax, alphaMin=.3, lineWidth=50, lineColor='#03045e')
lndF.plotMigrationNetwork(fig, ax, alphaMin=.3, lineWidth=35, lineColor='#03045e')
srv.plotClean(fig, ax, frame=False, bbox=bbox)
fig.savefig(
    path.join(OUT_PTH, '{}_{}_CLN.png'.format(LND_TYPE, ID)), 
    facecolor='w', bbox_inches='tight', pad_inches=0, dpi=300
)
###############################################################################
# GA Settings
############################################################################### 
(weightMale, weightFemale) = (.5, 1)
POP_SIZE = int(10*(lndM.trapsNumber*1.25))
(GENS, MAT, MUT, SEL, VERBOSE) = (
    500,
    {'mate': .3, 'cxpb': 0.5}, 
    {'mean': 0, 'sd': min([i[1]-i[0] for i in bbox])/5, 'mutpb': .4, 'ipb': .5},
    {'tSize': 3},
    True
)
lndM_GA = deepcopy(lndM)
lndF_GA = deepcopy(lndF)
# Needed auxiliary variables --------------------------------------------------
bbox = lndM.getBoundingBox()
trpMsk = srv.genFixedTrapsMask(lndM.trapsFixed)
###############################################################################
# Registering Functions for GA
############################################################################### 
((lndM_GA, lndF_GA), logbook) = srv.optimizeTwoSexesTrapsGA(
        lndM_GA, lndF_GA, sexWeights={'M': .5, 'F': .5}, 
        pop_size='auto', generations=GENS,
        mating_params=MAT, mutation_params=MUT, selection_params=SEL,
        fitFuns={'outer': np.mean, 'inner': np.max}, verbose=True
    )
srv.exportLog(logbook, OUT_PTH, '{}_{}_LOG'.format(LND_TYPE, ID))
###############################################################################
# Plot Landscape
############################################################################### 
(fig, ax) = plt.subplots(1, 1, figsize=(15, 15), sharey=False)
lndM.plotSites(fig, ax, size=100)
# Plot Networks ---------------------------------------------------------------
lndM_GA.plotMigrationNetwork(fig, ax, alphaMin=.3, lineWidth=50, lineColor='#03045e')
lndF_GA.plotMigrationNetwork(fig, ax, alphaMin=.3, lineWidth=35, lineColor='#03045e')
lndM_GA.plotTraps(fig, ax, colors={0: '#f7258522'}, lws=(2, 0), fill=True, ls='--', zorder=(25, 4))
lndF_GA.plotTraps(fig, ax, colors={0: '#a06cd522'}, lws=(2, 0), fill=True, ls=':', zorder=(25, 4))
srv.plotClean(fig, ax, frame=False, bbox=bbox)
srv.plotFitness(fig, ax, min(logbook['min']))
fig.savefig(
    path.join(OUT_PTH, '{}_{}_TRP.png'.format(LND_TYPE, ID)), 
    facecolor='w', bbox_inches='tight', pad_inches=0, dpi=300
)
plt.close('all')
print("Done!")