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


(OUT_PTH, LND_TYPE, ID, PTS_TYPE) = ('./scratch/', 'AGG', '001', 3)
ptsNum = 200
bbox = ((-150, 150), (-100, 100))
srv.makeFolder(OUT_PTH)
###############################################################################
# Point-process
###############################################################################
xy = srv.ptsRandUniform(ptsNum, bbox).T
pType = np.random.choice(PTS_TYPE, xy.shape[1])
points = pd.DataFrame({'x': xy[0], 'y': xy[1], 't': pType})
###############################################################################
# Mosquito Movement
###############################################################################
mKer = {'params': [.075, 1.0e-10, math.inf], 'zeroInflation': .75}
mMsk = np.array([
    [0.2, 0.8, 0.0],
    [0.0, 0.0, 1.0],
    [0.9, 0.1, 0.0]
])
###############################################################################
# Defining Traps
###############################################################################
traps = pd.DataFrame({
    'x': [-100, 50, -85, 75],
    'y': [75, -50, -75, 75],
    't': [1, 0, 1, 0],
    'f': [0, 0, 0, 0]
})
tKer = {
    0: {'kernel': srv.exponentialDecay, 'params': {'A': 0.50, 'b': .050}},
    1: {'kernel': srv.exponentialDecay, 'params': {'A': 0.35, 'b': .050}}
}
tMsk = np.asarray([
    [1, 0.0, 0.0],
    [0.0, 0.5, 0.5]
])
# tMsk = np.asarray([
#     [1.0, 1.0, 1.0],
#     [1.0, 1.0, 1.0]
# ])
###############################################################################
# Landscape
###############################################################################
lnd = srv.Landscape(
    points, 
    kernelParams=mKer, maskingMatrix=mMsk,
    traps=traps, trapsKernels=tKer, trapsMask=tMsk
)
###############################################################################
# Plotting
###############################################################################
bbox = lnd.getBoundingBox()
(fig, ax) = plt.subplots(1, 1, figsize=(15, 15), sharey=False)
lnd.plotSites(fig, ax, size=100)
lnd.plotMigrationNetwork(fig, ax, alphaMin=.6, lineWidth=25)
lnd.plotTraps(fig, ax)
lnd.plotTrapsNetwork(fig, ax)
srv.plotClean(fig, ax, frame=False)
fig.savefig(
    path.join(OUT_PTH, 'TRP_DEV.png'), 
    facecolor='w', bbox_inches='tight', pad_inches=0.0, dpi=300
)