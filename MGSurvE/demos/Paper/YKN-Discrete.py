#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
from os import path
from sys import argv
from copy import deepcopy
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import MGSurvE as srv


if srv.isNotebook():
    (ID, AP, RID) = ('YKN', 'man', '001')
else:
    (ID, AP, RID) = argv[1:]
PRINT_BLANK = False
###############################################################################
# File ID
###############################################################################
GENS = 10
OUT_PTH = './sims_out/'
srv.makeFolder(OUT_PTH)
###############################################################################
# File ID
###############################################################################
LND_PTH = './GEO/{}_LatLon.csv'.format(ID)
TRAP_TYP = [0]*8 + [1]*8
TRPS_NUM = len(TRAP_TYP)
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
    'kernelParams': {'params': srv.AEDES_EXP_PARAMS, 'zeroInflation': 1-0.28}
}
###############################################################################
# Defining Traps
###############################################################################
nullTraps = [0]*TRPS_NUM
cntr = ([np.mean(YK_LL['lon'])]*TRPS_NUM, [np.mean(YK_LL['lat'])]*TRPS_NUM)
sid = [0]*TRPS_NUM
traps = pd.DataFrame({
    'sid': sid,
    'lon': cntr[0], 'lat': cntr[1], 
    't': TRAP_TYP, 'f': nullTraps
})
# Setup trap kernels ----------------------------------------------------------
tKer = {
    1: {
        'kernel': srv.sigmoidDecay,     
        'params': {'A': 0.5, 'rate': .25, 'x0': 1/0.0629534}
    },
    0: {
        'kernel': srv.exponentialDecay, 
        'params': {'A': 0.5, 'b': 0.0629534}
    }
}
# meanDistances = [srv.nSolveKernel(tKer[i], 0.5, 20) for i in tKer.keys()]
###############################################################################
# Setting Landscape Up
###############################################################################
lnd = srv.Landscape(
    YK_LL, 
    kernelFunction=mKer['kernelFunction'], kernelParams=mKer['kernelParams'],
    traps=traps, trapsKernels=tKer, trapsRadii=[.50, .45, .40],
    landLimits=YK_BBOX
)
bbox = lnd.getBoundingBox()
lndGA = deepcopy(lnd)
###############################################################################
# Plot Landscape
###############################################################################
if PRINT_BLANK:
    (fig, ax) = (
        plt.figure(figsize=(15, 15)),
        plt.axes(projection=ccrs.PlateCarree())
    )
    lnd.plotSites(fig, ax, size=50)
    lnd.plotMigrationNetwork(
        fig, ax, 
        lineWidth=7.5, alphaMin=.05, alphaAmplitude=7.5
    )
    # lnd.plotLandBoundary(fig, ax)
    srv.plotClean(fig, ax, bbox=lnd.landLimits)
    fig.savefig(
        path.join(OUT_PTH, '{}D-{}_{:02d}-{:02d}_CLN.png'.format(ID, AP, TRPS_NUM, RID)), 
        facecolor='w', bbox_inches='tight', pad_inches=0.1, dpi=300
    )
    plt.close('all')
###############################################################################
# Registering Functions for GA
###############################################################################
if (AP=='man'):
    outer = np.mean
elif (AP=='sum'):
    outer = np.sum
elif (AP=='max'):
    outer = np.max
# Optimize discrete -----------------------------------------------------------
(lnd, logbook) = srv.optimizeDiscreteTrapsGA(
    lndGA, verbose=False,
    generations=GENS,
    pop_size='auto',
    mating_params='auto', 
    mutation_params='auto', 
    selection_params='auto',
    fitFuns={'inner': np.sum, 'outer': outer}
)
srv.exportLog(logbook, OUT_PTH, '{}D-{}_{:02d}-{:02d}_LOG'.format(ID, AP, TRPS_NUM, RID))
srv.dumpLandscape(lnd, OUT_PTH, '{}D-{}_{:02d}-{:02d}_TRP'.format(ID, AP, TRPS_NUM, RID), fExt='pkl')
###############################################################################
# Plots
###############################################################################
# Landscape -------------------------------------------------------------------
lnd = srv.loadLandscape(
    OUT_PTH, '{}D-{}_{:02d}-{:02d}_TRP'.format(ID, AP, TRPS_NUM, RID), 
    fExt='pkl'
)
(fig, ax) = (
    plt.figure(figsize=(15, 15)), 
    plt.axes(projection=ccrs.PlateCarree())
)
lnd.plotSites(fig, ax, size=50)
# lnd.plotMigrationNetwork(fig, ax, lineWidth=7.5, alphaMin=.05, alphaAmplitude=7.5)
lnd.plotTraps(fig, ax, zorders=(30, 25))
# srv.plotFitness(fig, ax, min(logbook['min']), fmt='{:.5f}', fontSize=100)
srv.plotClean(fig, ax, bbox=lnd.landLimits)
fig.savefig(
    path.join(OUT_PTH, '{}D-{}_{:02d}-{:02d}_TRP.png'.format(ID, AP, TRPS_NUM, RID)), 
    facecolor='w', bbox_inches='tight', pad_inches=0.1, dpi=300
)
plt.close('all')
# Traps Kernels ---------------------------------------------------------------
(fig, ax) = plt.subplots(1, 1, figsize=(15, 5), sharey=False)
(fig, ax) = srv.plotTrapsKernels(fig, ax, lnd, distRange=(0, 100), aspect=.175)
fig.savefig(
    path.join(OUT_PTH, '{}D-{}_{:02d}-{:02d}_KER.png'.format(ID, AP, TRPS_NUM, RID)), 
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
    path.join(OUT_PTH, '{}D-{}_{:02d}-{:02d}_GA.png'.format(ID, AP, TRPS_NUM, RID)), 
    facecolor='w', bbox_inches='tight', pad_inches=0.1, dpi=300
)
plt.close('all')



