#!/usr/bin/env python
# -*- coding: utf-8 -*-

import warnings
import numpy as np
import pandas as pd
from os import path
from sys import argv
from copy import deepcopy
import matplotlib.pyplot as plt
from compress_pickle import dump, load
import Constants as cst
import MGSurvE as srv
warnings.filterwarnings('ignore', 'The iteration is not making good progress')


(GENS, VERBOSE, OUT_PTH) = (cst.gens, cst.verbose, cst.out_pth)
if srv.isNotebook():
    (ID, ZIK, COMBO) = ('Grid_LND_HOM', 'ZI', 'sum-men')
else:
    (ID, ZIK, COMBO) = (argv[1], argv[2], argv[3])
###############################################################################
# Params for Combos
###############################################################################
PREP = 'LngHi-'
fCombo = {
    'sum-men': [np.sum,  np.mean],
    'men-max': [np.mean, np.max ],
    'max-men': [np.max,  np.mean],
    'max-max': [np.max,  np.max ],
    'men-men': [np.mean, np.mean],
    'men-sum': [np.mean, np.sum ],
    'sum-max': [np.sum,  np.max ],
    'max-sum': [np.max,  np.sum ],
}
###############################################################################
# Load Landscape
###############################################################################
ID = f'{ZIK}-{ID}'
lnd = srv.loadLandscape(OUT_PTH, ID)
# Needed auxiliary variables --------------------------------------------------
(bbox, trpMsk) = (lnd.getBoundingBox(), srv.genFixedTrapsMask(lnd.trapsFixed))
###############################################################################
# GA Settings
############################################################################### 
TRPS_NUM = lnd.trapsCoords.shape[0]
POP_SIZE = int(10*(lnd.trapsNumber*1.5))
(MAT, MUT, SEL) = cst.gaParams
lndGA = deepcopy(lnd)
###############################################################################
# Registering Functions for GA
###############################################################################
fitFuns = {'inner': fCombo[COMBO][0], 'outer': fCombo[COMBO][1]}
(lnd, logbook) = srv.optimizeDiscreteTrapsGA(
    lndGA, pop_size=POP_SIZE, generations=GENS,
    mating_params=MAT, mutation_params=MUT, selection_params=SEL,
    fitFuns=fitFuns, verbose=VERBOSE
)
srv.dumpLandscape(lnd, OUT_PTH, PREP+'{}_TRP-{}'.format(ID, COMBO), fExt='pkl')
srv.exportLog(logbook, OUT_PTH, PREP+'{}_TRP-{}'.format(ID, COMBO))
###############################################################################
# Plot GA
############################################################################### 
(fig, ax) = plt.subplots(figsize=(15, 15))
(fig, ax) = srv.plotGAEvolution(fig, ax, logbook)
pthSave = path.join(OUT_PTH, PREP+'{}_GAP-{}'.format(ID, COMBO))
fig.savefig(
    pthSave,
    facecolor='w', bbox_inches='tight', 
    pad_inches=.1, dpi=300
)
# Export plots ----------------------------------------------------------------
lnd = srv.loadLandscape(OUT_PTH, PREP+'{}_TRP-{}'.format(ID, COMBO), fExt='pkl')
bbox = lnd.getBoundingBox()
trpMsk = srv.genFixedTrapsMask(lnd.trapsFixed)
(fig, ax) = plt.subplots(1, 1, figsize=(15, 15), sharey=False)
lnd.plotSites(fig, ax, size=200)
lnd.plotMaskedMigrationNetwork(fig, ax, alphaMin=.5, lineWidth=50)
lnd.plotTraps(fig, ax, size=200)
srv.plotClean(fig, ax, bbox=bbox, frame=False, pad=cst.pad_i)
# srv.plotFitness(fig, ax, min(logbook['min']), zorder=30)
fig.savefig(
    path.join(OUT_PTH, PREP+'{}_TRP-{}.png'.format(ID, COMBO)), 
    facecolor='w', bbox_inches='tight', pad_inches=cst.pad, dpi=cst.dpi
)
plt.close('all')