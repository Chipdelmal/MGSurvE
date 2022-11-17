#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
from os import path
from copy import deepcopy
import cartopy.crs as crs
import matplotlib.pyplot as plt
import MGSurvE as srv

ID = 'TTP'
###############################################################################
# File ID
###############################################################################
GENS = 1500
OUT_PTH = './sims_out/'
srv.makeFolder(OUT_PTH)
###############################################################################
# File ID
###############################################################################
LND_PTH = './GEO/{}_LatLon.csv'.format(ID)
TRPS_NUM = 10
TRAP_TYP = [0, 0, 1, 0, 1, 1, 0, 1, 0, 2]
###############################################################################
# Load pointset
###############################################################################
YK_LL = pd.read_csv(LND_PTH, names=['lon', 'lat'])
YK_LL['t'] = [0]*YK_LL.shape[0]
pad = 0.00125
YK_BBOX = (
    (min(YK_LL['lon'])-pad, max(YK_LL['lon'])+pad),
    (min(YK_LL['lat'])-pad, max(YK_LL['lat'])+pad)
)
# YK_LL = YK_LL.reindex(columns=['lat', 'lon'])
# Movement Kernel -------------------------------------------------------------
mKer = {
    'kernelFunction': srv.zeroInflatedExponentialKernel,
    'kernelParams': {'params': srv.AEDES_EXP_PARAMS, 'zeroInflation': 0}
}
###############################################################################
# Defining Traps
###############################################################################
nullTraps = [0]*TRPS_NUM
cntr = ([np.mean(YK_LL['lon'])]*TRPS_NUM, [np.mean(YK_LL['lat'])]*TRPS_NUM)
traps = pd.DataFrame({
    'lon': cntr[0], 'lat': cntr[1], 
    't': TRAP_TYP, 'f': nullTraps
})
# Setup trap kernels ----------------------------------------------------------
tKer = {
    2: {
        'kernel': srv.sigmoidDecay,     
        'params': {'A': .8, 'rate': .15, 'x0': 30}
    },
    1: {
        'kernel': srv.exponentialDecay, 
        'params': {'A': 1, 'b': 0.0425}
    },
    0: {
        'kernel': srv.exponentialAttractiveness,
        'params': {'A': 1, 'k': .01, 's': .3, 'gamma': .975, 'epsilon': 0}
    }
}
###############################################################################
# Setting Landscape Up
###############################################################################
lnd = srv.Landscape(
    YK_LL, 
    kernelFunction=mKer['kernelFunction'], kernelParams=mKer['kernelParams'],
    traps=traps, trapsKernels=tKer, trapsRadii=[.9, .8, .75],
    landLimits=YK_BBOX
)
bbox = lnd.getBoundingBox()
lndGA = deepcopy(lnd)
###############################################################################
# Registering Functions for GA
############################################################################### 
(lnd, logbook) = srv.optimizeTrapsGA(
    lndGA, pop_size='auto', generations=GENS,
    mating_params='auto', mutation_params='auto', selection_params='auto',
    fitFuns={'outer': np.mean, 'inner': np.max}
)
srv.exportLog(logbook, OUT_PTH, '{}_LOG'.format(ID))
srv.dumpLandscape(lnd, OUT_PTH, '{}_{:02d}_TRP'.format(ID, TRPS_NUM), fExt='pkl')
###############################################################################
# Plots
###############################################################################
# Landscape -------------------------------------------------------------------
lnd = srv.loadLandscape(OUT_PTH, '{}_{:02d}_TRP'.format(ID, TRPS_NUM), fExt='pkl')
(fig, ax) = (
    plt.figure(figsize=(15, 15)), plt.axes(projection=crs.PlateCarree())
)
lnd.plotSites(fig, ax, size=50)
lnd.plotMigrationNetwork(fig, ax, lineWidth=20, alphaMin=.075, alphaAmplitude=7.5)
lnd.plotTraps(fig, ax, zorders=(30, 25))
# srv.plotFitness(fig, ax, min(logbook['min']), fmt='{:.5f}', fontSize=100)
srv.plotClean(fig, ax, bbox=YK_BBOX)
fig.savefig(
    path.join(OUT_PTH, '{}_{:02d}_TRP.png'.format(ID, TRPS_NUM)), 
    facecolor='w', bbox_inches='tight', pad_inches=0.1, dpi=300
)
plt.close('all')
# Traps Kernels ---------------------------------------------------------------
(fig, ax) = plt.subplots(1, 1, figsize=(15, 5), sharey=False)
(fig, ax) = srv.plotTrapsKernels(fig, ax, lnd, distRange=(0, 100), aspect=.175)
fig.savefig(
    path.join(OUT_PTH, '{}_{:02d}_KER.png'.format(ID, TRPS_NUM)), 
    facecolor='w', bbox_inches='tight', pad_inches=0.1, dpi=300
)
plt.close('all')
# GA --------------------------------------------------------------------------
log = pd.DataFrame(logbook)
log.rename(columns={'median': 'avg'}, inplace=True)
(fig, ax) = plt.subplots(1, 1, figsize=(15, 5), sharey=False)
srv.plotGAEvolution(
    fig, ax,
    logbook,
    colors={'mean': '#ffffff', 'envelope': '#1565c0'},
    alphas={'mean': .1, 'envelope': 0.5},
    aspect=1/3
)
ax.set_ylim(-10, 5000)
ax.set_aspect((1/3)/ax.get_data_ratio())
fig.savefig(
     path.join(OUT_PTH, '{}_{:02d}_GA.png'.format(ID, TRPS_NUM)), 
    facecolor='w', bbox_inches='tight', pad_inches=0.1, dpi=300
)
plt.close('all')


