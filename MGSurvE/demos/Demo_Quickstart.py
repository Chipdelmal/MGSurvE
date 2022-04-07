#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
from os import path
from copy import deepcopy
import matplotlib.pyplot as plt
import MGSurvE as srv


(ID, OUT_PTH) = ('QSTART', './QSTART/')
srv.makeFolder(OUT_PTH)
###############################################################################
# Defining Landscape (Random in a Ring)
###############################################################################
ptsNum = 150
radii = (75, 100)
xy = srv.ptsDonut(ptsNum, radii).T
points = pd.DataFrame({'x': xy[0], 'y': xy[1], 't': [0]*ptsNum})
###############################################################################
# Defining Traps
###############################################################################
nullTraps = [0, 0, 0, 0, 0]
traps = pd.DataFrame({
    'x': nullTraps, 'y': nullTraps,
    't': nullTraps, 'f': nullTraps
})
tKer = {0: {'kernel': srv.exponentialDecay, 'params': {'A': .5, 'b': .1}}}
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
lnd.plotMigrationNetwork(fig, ax, alphaMin=.6, lineWidth=25)
lnd.plotTraps(fig, ax)
srv.plotClean(fig, ax, frame=False)
srv.plotFitness(fig, ax, min(logbook['min']))
fig.savefig(
    path.join(OUT_PTH, '{}_TRP.png'.format(ID)), 
    facecolor='w', bbox_inches='tight', pad_inches=0.1, dpi=300
)
plt.close('all')
