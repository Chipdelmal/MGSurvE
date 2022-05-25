#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
from os import path
from copy import deepcopy
import matplotlib.pyplot as plt
import MGSurvE as srv


(ID, OUT_PTH) = ('qstart', './demos_out/')
srv.makeFolder(OUT_PTH)
###############################################################################
# Defining Landscape (Random in a Ring)
###############################################################################
ptsNum = 150
radii = (75, 100)
xy = srv.ptsDonut(ptsNum, radii).T
pType = np.random.choice(1, xy.shape[1])
points = pd.DataFrame({'x': xy[0], 'y': xy[1], 't': pType})
###############################################################################
# Defining Traps
###############################################################################
nullTraps = [0, 0, 0, 0, 0]
traps = pd.DataFrame({
    'x': nullTraps, 'y': nullTraps,
    't': [0, 0, 0, 1, 1], 'f': nullTraps
})
tKer = {
    0: {'kernel': srv.exponentialDecay, 'params': {'A': .5, 'b': .1}},
    1: {'kernel': srv.exponentialDecay, 'params': {'A': .5, 'b': .125}}
}
###############################################################################
# Setting Landscape Object Up
###############################################################################
lnd = srv.Landscape(
    points, 
    kernelParams={'params': srv.MEDIUM_MOV_EXP_PARAMS, 'zeroInflation': .25}, 
    traps=traps, trapsKernels=tKer
)
###############################################################################
# Plotting Traps' Positions
############################################################################### 
lndGA = deepcopy(lnd)
(lnd, logbook) = srv.optimizeTrapsGA(
    lndGA, generations=500, 
    pop_size='auto', mating_params='auto', 
    mutation_params='auto', selection_params='auto',
)
srv.exportLog(logbook, OUT_PTH, '{}_LOG'.format(ID))
###############################################################################
# Plotting Optimized Landscape
############################################################################### 
(fig, ax) = plt.subplots(1, 1, figsize=(15, 15), sharey=False)
lnd.plotSites(fig, ax, size=100)
lnd.plotMaskedMigrationNetwork(fig, ax, alphaMin=.6, lineWidth=25)
lnd.plotTraps(fig, ax)# , zorders=(60, 50))
srv.plotClean(fig, ax, frame=False)
srv.plotFitness(fig, ax, min(logbook['min']))
fig.savefig(
    path.join(OUT_PTH, '{}_TRP.png'.format(ID)), 
    facecolor='w', bbox_inches='tight', pad_inches=0.1, dpi=300
)
plt.close('all')
###############################################################################
# Plotting Optimized Landscape
############################################################################### 
(fig, ax) = plt.subplots(1, 1, figsize=(15, 15), sharey=False)
srv.plotMatrix(fig, ax, lnd.trapsMigration, vmax=1e-2, trapsNumber=len(nullTraps))
srv.plotClean(fig, ax, frame=False)
fig.savefig(
    path.join(OUT_PTH, '{}_MTX.png'.format(ID)), 
    facecolor='w', bbox_inches='tight', pad_inches=0.1, dpi=300
)
plt.close('all')