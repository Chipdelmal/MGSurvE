#!/usr/bin/env python
# -*- coding: utf-8 -*-

from sys import argv
import math
import numpy as np
import pandas as pd
from os import path
from sys import argv
from copy import deepcopy
import matplotlib.pyplot as plt
from deap import base, creator, algorithms, tools
from compress_pickle import dump, load
from sklearn.preprocessing import normalize
import MGSurvE as srv
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import warnings
warnings.filterwarnings('ignore', 'The iteration is not making good progress')

(FXD_TRPS, TRPS_NUM, GENS) = (False, 10, 100)
DIAG_VAL = 0.01
###############################################################################
# Debugging fixed traps at land masses
###############################################################################
OUT_PTH = './sims_out/'
if FXD_TRPS:
    ID = 'STP_FXD'
else:
    ID = 'STP_FXN'
###############################################################################
# Load Pointset
###############################################################################
SAO_TOME_LL = pd.read_csv(path.join('./GEO', 'STP_LatLon.csv'))
SAO_TOME_LL['t'] = [0]*SAO_TOME_LL.shape[0]
SAO_bbox = (
    (min(SAO_TOME_LL['x']), max(SAO_TOME_LL['x'])),
    (min(SAO_TOME_LL['y']), max(SAO_TOME_LL['y']))
)
SAO_cntr = [i[0]+(i[1]-i[0])/2 for i in SAO_bbox]
SAO_LIMITS = ((6.41, 6.79), (-0.0475, .45))
# Get location of minor land-masses -------------------------------------------
SAO_FIXED = [tuple(SAO_TOME_LL.loc[i][['x', 'y']]) for i in (24, 212)]
FXD_NUM = len(SAO_FIXED)
###############################################################################
# Load Migration Matrix
###############################################################################
migration = np.genfromtxt(
    path.join('./GEO', 'STP_Migration.csv'), delimiter=','
)
np.fill_diagonal(migration, DIAG_VAL)
SAO_TOME_MIG = normalize(migration, axis=1, norm='l1')
###############################################################################
# Defining Traps
#   North and South land masses (mini islands): 51, 239 (zero-indexed)
###############################################################################
(initTyp, initFxd) = ([0]*TRPS_NUM, [0]*TRPS_NUM)
(initLon, initLat) = ([SAO_cntr[0]]*TRPS_NUM, [SAO_cntr[1]]*TRPS_NUM)
if FXD_TRPS:
    for i in range(FXD_NUM):
        initFxd[TRPS_NUM-(i+1)] = 1
        initLon[TRPS_NUM-(i+1)] = SAO_FIXED[i][0]
        initLat[TRPS_NUM-(i+1)] = SAO_FIXED[i][1]
traps = pd.DataFrame({'x': initLon, 'y': initLat, 't': initTyp, 'f': initFxd})
tKer = {0: {'kernel': srv.exponentialDecay, 'params': {'A': .5, 'b': 100}}}
###############################################################################
# Setting Landscape Up
###############################################################################
lnd = srv.Landscape(
    SAO_TOME_LL, migrationMatrix=SAO_TOME_MIG,
    traps=traps, trapsKernels=tKer, landLimits=SAO_LIMITS
)
bbox = lnd.getBoundingBox()
trpMsk = srv.genFixedTrapsMask(lnd.trapsFixed)
###############################################################################
# Plot Results
###############################################################################
(fig, ax) = (
    plt.figure(figsize=(15, 15)),
    plt.axes(projection=ccrs.PlateCarree())
)
lnd.plotSites(fig, ax, size=250)
lnd.plotMigrationNetwork(
    fig, ax, lineWidth=10, alphaMin=.1, alphaAmplitude=2.5
)
lnd.plotLandBoundary(fig, ax)
srv.plotClean(fig, ax, bbox=lnd.landLimits)
fig.savefig(
    path.join(OUT_PTH, '{}_{:02d}_CLN.png'.format(ID, TRPS_NUM)), 
    facecolor='w', bbox_inches='tight', pad_inches=0.1, dpi=400
)
plt.close('all')
###############################################################################
# GA Settings
############################################################################### 
POP_SIZE = int(10*(lnd.trapsNumber*1.25))
(MAT, MUT, SEL) = (
    {'mate': .35, 'cxpb': 0.5}, 
    {
        'mean': 0, 'sd': max([abs(i[1]-i[0]) for i in bbox])/5, 
        'mutpb': .35, 'ipb': .5
    },
    {'tSize': 5}
)
lndGA = deepcopy(lnd)
# Reducing the bbox for init sampling -----------------------------------------
redFract = 0
reduction = [(i[1]-i[0])/2*redFract for i in bbox]
bboxRed = [(i[0]+r, i[1]-r) for (i, r) in zip(bbox,reduction)]
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
# Plot Results
###############################################################################
(fig, ax) = (
    plt.figure(figsize=(15, 15)),
    plt.axes(projection=ccrs.PlateCarree())
)
lnd.plotSites(fig, ax, size=250)
lnd.plotMigrationNetwork(
    fig, ax, lineWidth=10, alphaMin=.1, alphaAmplitude=2.5
)
lnd.plotTraps(fig, ax, zorders=(25, 20))
srv.plotFitness(fig, ax, min(logbook['min']), fmt='{:.2f}')
lnd.plotLandBoundary(fig, ax)
srv.plotClean(fig, ax, bbox=lnd.landLimits)
fig.savefig(
    path.join(OUT_PTH, '{}_{:02d}_TRP.png'.format(ID, TRPS_NUM)), 
    facecolor='w', bbox_inches='tight', pad_inches=0.1, dpi=400
)
plt.close('all')


