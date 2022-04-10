#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math
import numpy as np
import pandas as pd
from os import path
from copy import deepcopy
import matplotlib.pyplot as plt
import MGSurvE as srv


(ID, OUT_PTH) = ('GA_DEMO_PTS', './demos_out/')
(PTS_NUM, PTS_TYPE, LND) = (150, 3, 'HOM')
srv.makeFolder(OUT_PTH)
###############################################################################
# Generate Sites
###############################################################################
bside = math.sqrt(PTS_NUM)*11.25
bbox = ((-bside, bside), (-bside, bside))
xy = srv.ptsRegularGrid(int(math.sqrt(PTS_NUM)), bbox).T
if LND == 'HOM':
    pType = [0]*xy.shape[1]
else:
    pType = np.random.choice(PTS_TYPE, xy.shape[1])
points = pd.DataFrame({'x': xy[0], 'y': xy[1], 't': pType})
###############################################################################
# Setup Points Movement
###############################################################################
msk = [
    [0.100, 0.700, 0.200],
    [0.100, 0.100, 0.800],
    [0.750, 0.125, 0.125],
]
movKer = {'params': [.075, 1.0e-10, math.inf], 'zeroInflation': .75}
###############################################################################
# Generate Traps
###############################################################################
trp = [
    [0, 0, 0, 0],
    [0, 0, 0, 0],
    [0, 0, 0, 0],
]
traps = pd.DataFrame(trp, columns=['x', 'y', 't', 'f'])
tker = {0: {'kernel': srv.exponentialDecay, 'params': {'A': 0.5, 'b': .05}}}
###############################################################################
# Setup Landscapes
###############################################################################
if LND == 'HOM':
    lnd = srv.Landscape(
        points, traps=traps, trapsKernels=tker
    )
else:
    lnd = srv.Landscape(
        points, maskingMatrix=msk, traps=traps, trapsKernels=tker
    )
bbox = lnd.getBoundingBox()
trpMsk = srv.genFixedTrapsMask(lnd.trapsFixed)
###############################################################################
# Plot Landscapes
###############################################################################
(fig, ax) = plt.subplots(1, 1, figsize=(10, 10), sharey=False)
lnd.plotSites(fig, ax, size=100)
lnd.plotMaskedMigrationNetwork(fig, ax, alphaMin=.6, lineWidth=25)
srv.plotClean(fig, ax, frame=False)
fig.savefig(
    path.join(OUT_PTH, '{}_{}_CLN.png'.format(ID, LND)), 
    facecolor='w', bbox_inches='tight', pad_inches=0.1, dpi=300
)
plt.close('all')
###############################################################################
# GA
###############################################################################
POP_SIZE = int(10*(lnd.trapsNumber*1.25))
(GENS, MAT, MUT, SEL) = (
    500,
    {'mate': .3, 'cxpb': 0.5}, 
    {'mean': 0, 'sd': min([i[1]-i[0] for i in bbox])/5, 'mutpb': .4, 'ipb': .5},
    {'tSize': 3}
)
lndGA = deepcopy(lnd)
###############################################################################
# Registering Functions for GA
############################################################################### 
(lnd, logbook) = srv.optimizeTrapsGA(
        lndGA, pop_size='auto', generations=GENS,
        mating_params=MAT, mutation_params=MUT, selection_params=SEL,
        fitFuns={'outer': np.mean, 'inner': np.max}, verbose=True
    )
srv.exportLog(logbook, OUT_PTH, '{}_LOG'.format(ID))
###############################################################################
# Plot Landscape
############################################################################### 
(fig, ax) = plt.subplots(1, 1, figsize=(15, 15), sharey=False)
lnd.plotSites(fig, ax, size=100)
lnd.plotMaskedMigrationNetwork(fig, ax, alphaMin=.6, lineWidth=25)
lnd.plotTraps(fig, ax)
srv.plotClean(fig, ax, frame=False)
srv.plotFitness(fig, ax, min(logbook['min']))
fig.savefig(
    path.join(OUT_PTH, '{}_{}_TRP'.format(ID, LND)), 
    facecolor='w', bbox_inches='tight', pad_inches=0.1, dpi=300
)
plt.close('all')