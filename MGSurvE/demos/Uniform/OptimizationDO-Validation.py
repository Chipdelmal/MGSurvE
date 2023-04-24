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

import os
os.environ["OMP_NUM_THREADS"] = "24" # export OMP_NUM_THREADS=4
os.environ["OPENBLAS_NUM_THREADS"] = "24" # export OPENBLAS_NUM_THREADS=4 
os.environ["MKL_NUM_THREADS"] = "24" # export MKL_NUM_THREADS=6
os.environ["VECLIB_MAXIMUM_THREADS"] = "24" # export VECLIB_MAXIMUM_THREADS=4
os.environ["NUMEXPR_NUM_THREADS"] = "24" # export NUMEXPR_NUM_THREADS=6


(GENS, VERBOSE, OUT_PTH) = (cst.gens, cst.verbose, cst.out_pth)
if srv.isNotebook():
    (ID, ZIK, COMBO, RID) = ('Grid_LND_HOM', 'ZI', 'sum-men', '1')
else:
    (ID, ZIK, COMBO, RID) = (argv[1], argv[2], argv[3], argv[4])
RID = int(RID)
###############################################################################
# Params for Combos
###############################################################################
PREP = 'ShtLo-'
fCombo = {
    'sum-men': [np.sum,  np.mean],
    'sum-max': [np.sum,  np.max ],
    'sum-sum': [np.sum,  np.sum ],
    'sum-min': [np.sum,  np.min ],
    'men-men': [np.mean, np.mean],
    'men-max': [np.mean, np.max ],
    'men-sum': [np.mean, np.sum ],
    'men-min': [np.mean, np.min ],
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
# print("*"*100)
# print("* Running {}".format(COMBO))
# print("*"*100)
fitFuns = {'inner': fCombo[COMBO][0], 'outer': fCombo[COMBO][1]}
(lnd, logbook) = srv.optimizeDiscreteTrapsGA(
    lndGA, pop_size=POP_SIZE, generations=GENS,
    mating_params=MAT, mutation_params=MUT, selection_params=SEL,
    fitFuns=fitFuns, verbose=VERBOSE
)
srv.dumpLandscape(lnd, OUT_PTH, PREP+'{}_TRP-{}-{:02d}'.format(ID, COMBO, RID), fExt='pkl')
srv.exportLog(logbook, OUT_PTH, PREP+'{}_TRP-{}-{:02d}'.format(ID, COMBO, RID))
###############################################################################
# Plot GA
############################################################################### 
(fig, ax) = plt.subplots(figsize=(15, 15))
(fig, ax) = srv.plotGAEvolution(fig, ax, logbook)
pthSave = path.join(OUT_PTH, PREP+'{}_GAP-{}-{:02d}'.format(ID, COMBO, RID))
fig.savefig(
    pthSave,
    facecolor='w', bbox_inches='tight', 
    pad_inches=.1, dpi=300
)
# Export plots ----------------------------------------------------------------
lnd = srv.loadLandscape(OUT_PTH, PREP+'{}_TRP-{}-{:02d}'.format(ID, COMBO, RID), fExt='pkl')
bbox = lnd.getBoundingBox()
trpMsk = srv.genFixedTrapsMask(lnd.trapsFixed)
(fig, ax) = plt.subplots(1, 1, figsize=(15, 15), sharey=False)
lnd.plotSites(fig, ax, size=200)
lnd.plotMaskedMigrationNetwork(fig, ax, alphaMin=.5, lineWidth=50)
lnd.plotTraps(fig, ax, size=200)
srv.plotClean(fig, ax, bbox=bbox, frame=False, pad=cst.pad_i)
# srv.plotFitness(fig, ax, min(logbook['min']), zorder=30)
fig.savefig(
    path.join(OUT_PTH, PREP+'{}_TRP-{}-{:02d}'.format(ID, COMBO, RID)), 
    facecolor='w', bbox_inches='tight', pad_inches=cst.pad, dpi=cst.dpi
)
plt.close('all')