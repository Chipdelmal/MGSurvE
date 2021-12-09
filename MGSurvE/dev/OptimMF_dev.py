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

if srv.isNotebook():
    (OUT_PTH, LND_TYPE, ID) = ('./Lands', 'DNUT', 'D01')
else:
    (OUT_PTH, LND_TYPE, ID) = (argv[1], argv[2], argv[3].zfill(3))
###############################################################################
# Generating Pointsets
###############################################################################
if LND_TYPE == 'UNIF':
    ptsNum = 400
    bbox = ((-225, 225), (-175, 175))
    xy = srv.ptsRandUniform(ptsNum, bbox).T
elif LND_TYPE == 'GRID':
    ptsNum = 20
    bbox = ((-225, 225), (-225, 225))
    xy = srv.ptsRegularGrid(ptsNum, bbox).T
elif LND_TYPE == 'DNUT':
    ptsNum = 300
    radii = (100, 150)
    xy = srv.ptsDonut(ptsNum, radii).T
points = pd.DataFrame({'x': xy[0], 'y': xy[1], 't': [0]*xy.shape[1]})
###############################################################################
# Defining Traps
###############################################################################
traps = pd.DataFrame({
    'x': [0, 0, 0, 0, 0],
    'y': [0, 0, 0, 0, 0],
    't': [0, 0, 0, 0, 0],
    'f': [0, 0, 0, 0, 0]
})
tKernels = {
    'Male': {
        0: {'kernel': srv.exponentialDecay, 'params': {'A': .3, 'b': .05}}
    },
    'Female': {
        0: {'kernel': srv.exponentialDecay, 'params': {'A': .25, 'b': .025}}
    }
}
###############################################################################
# Defining Movement
###############################################################################
movementKernel = {
    'Male': {
        'kernelFunction': srv.zeroInflatedExponentialKernel,
        'kernelParams': {'params': [.075, 1.0e-10, math.inf], 'zeroInflation': .5}
    },
    'Female': {
        'kernelFunction': srv.zeroInflatedExponentialKernel,
        'kernelParams': {'params': [.075, 1.0e-10, math.inf], 'zeroInflation': .75}
    }
}
###############################################################################
# Setting Landscape Up
###############################################################################
lndMale = srv.Landscape(
    points, traps=traps,
    kernelFunction=movementKernel['Male']['kernelFunction'],
    kernelParams=movementKernel['Male']['kernelParams'],
    trapsKernels=tKernels['Male']
)
lndFemale = srv.Landscape(
    points, traps=traps,
    kernelFunction=movementKernel['Female']['kernelFunction'],
    kernelParams=movementKernel['Female']['kernelParams'],
    trapsKernels=tKernels['Female']
)