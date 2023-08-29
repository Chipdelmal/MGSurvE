#!/usr/bin/env python
# -*- coding: utf-8 -*-

CORES = 8
###############################################################################
# Load libraries and limit cores
###############################################################################
import os
os.environ["OMP_NUM_THREADS"] = str(CORES)
os.environ["OPENBLAS_NUM_THREADS"] = str(CORES)
os.environ["MKL_NUM_THREADS"] = str(CORES)
os.environ["VECLIB_MAXIMUM_THREADS"] = str(CORES)
os.environ["NUMEXPR_NUM_THREADS"] = str(CORES)
# Load libraries --------------------------------------------------------------
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
    (ID, ZIK) = ('Ring_LND_HOM', 'ZI')
else:
    (ID, ZIK) = (argv[1], argv[2])
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
(lnd, logbook) = srv.optimizeTrapsGA(
    lndGA, pop_size='auto', generations=GENS,
    mating_params=MAT, mutation_params=MUT, selection_params=SEL,
    fitFuns={'outer': np.max, 'inner': np.sum}, verbose=VERBOSE
)
srv.dumpLandscape(lnd, OUT_PTH, '{}_TRP-COS'.format(ID), fExt='pkl')
srv.exportLog(logbook, OUT_PTH, '{}_LOG-COS'.format(ID))
###############################################################################
# Plot GA
###############################################################################
logbook = pd.read_csv(path.join(OUT_PTH, '{}_LOG-COS.csv'.format(ID)))
log = pd.DataFrame(logbook)
(fig, ax) = plt.subplots(figsize=(15, 15))
(fig, ax) = srv.plotGAEvolution(fig, ax, logbook)
pthSave = path.join( OUT_PTH, '{}_GAP-COS'.format(ID))
fig.savefig(
    pthSave,
    facecolor='w', bbox_inches='tight', 
    pad_inches=.1, dpi=300
)
# Export plots ----------------------------------------------------------------
lnd = srv.loadLandscape(OUT_PTH, '{}_TRP-COS'.format(ID), fExt='pkl')
bbox = lnd.getBoundingBox()
(fig, ax) = plt.subplots(1, 1, figsize=(15, 15), sharey=False)
lnd.plotSites(fig, ax, size=200)
lnd.plotMaskedMigrationNetwork(fig, ax, alphaMin=.5, lineWidth=50)
lnd.plotTraps(fig, ax, size=200)
srv.plotClean(fig, ax, bbox=bbox, frame=False, pad=cst.pad_i)
# srv.plotFitness(fig, ax, min(logbook['min']), zorder=30)
fig.savefig(
    path.join(OUT_PTH, '{}_TRP-COS.png'.format(ID)), 
    facecolor='w', bbox_inches='tight', pad_inches=cst.pad, dpi=cst.dpi
)
plt.close('all')