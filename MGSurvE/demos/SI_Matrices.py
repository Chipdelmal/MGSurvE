#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math
from os import path
import pandas as pd
import MGSurvE as srv
import matplotlib.pyplot as plt


OUT_PTH = './demos_out/'
srv.makeFolder(OUT_PTH)
###############################################################################
# XY Landscape with two point types
###############################################################################
pts = [
    [0, 10,  0], 
    [-10, 0, 0], 
    [0, 20,  1],
    [-20, 30, 0]
]
points = pd.DataFrame(pts, columns=['x', 'y', 't'])
mskMigration = [
    [0.3, 0.7],
    [0.7, 0.3]
]
###############################################################################
# XY Landscape with two trap types
###############################################################################
trp = [
    [-2,  15, 1, 0],
    [-15, 30, 0, 0]
]
traps = pd.DataFrame(trp, columns=['x', 'y', 't', 'f'])
tker = {
    0: {'kernel': srv.exponentialDecay, 'params': {'A': 1, 'b': 0.5}},
    1: {'kernel': srv.sigmoidDecay,     'params': {'A': 1, 'rate': 0.75, 'x0': 2}}
}
###############################################################################
# Landscape Object
###############################################################################
lnd = srv.Landscape(
    points,
    distanceFunction=math.dist,
    kernelFunction=srv.zeroInflatedExponentialKernel,
    kernelParams={'params': srv.AEDES_EXP_PARAMS, 'zeroInflation': .72},
    maskingMatrix=mskMigration,
    traps=traps, trapsKernels=tker
)
lnd.distanceMatrix
lnd.migrationMatrix
lnd.maskedMigration
# Plot landscape --------------------------------------------------------------
(fig, ax) = plt.subplots(1, 1, figsize=(15, 15), sharey=False)
lnd.plotSites(fig, ax)
lnd.plotMaskedMigrationNetwork(fig, ax)
lnd.plotTraps(fig, ax)
# lnd.plotTrapsNetwork(fig, ax)
srv.plotClean(fig, ax, frame=True)
###############################################################################
# Plot
###############################################################################
(fig, ax) = plt.subplots(1, 2, figsize=(15, 15), sharey=False)
lnd.plotSites(fig, ax[0])
lnd.plotMaskedMigrationNetwork(fig, ax[0])
lnd.plotTraps(fig, ax[0])
lnd.plotTrapsNetwork(fig, ax[0])
srv.plotMatrix(fig, ax[1], lnd.trapsMigration, lnd.trapsNumber)
[srv.plotClean(fig, i, frame=False) for i in ax]