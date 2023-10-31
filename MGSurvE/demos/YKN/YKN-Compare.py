#!/usr/bin/env python
# -*- coding: utf-8 -*-

from glob import glob
import numpy as np
import pandas as pd
from os import path
from sys import argv
import cartopy.crs as ccrs
import shapely.geometry as sgeom
from cartopy.geodesic import Geodesic
import matplotlib.pyplot as plt
import MGSurvE as srv
import auxiliary as aux
from copy import deepcopy
plt.rcParams['axes.facecolor']='#00000000'
plt.rcParams['savefig.facecolor']='#00000000'


###############################################################################
# Bash and user inputs
###############################################################################
if srv.isNotebook():
    # User input (interactive session)
    (ID, AP, RID) = ('YKND', 'man', '01')
else:
    # Bash call input
    (ID, AP, RID) = argv[1:]
RID = int(RID)
FPAT = ID+'-{}_16*'
(MPATS, COLS, SCAL) =  (
    ['man', 'max'],
    ['#072ac877', '#8338ec77', '#f7258577'],
    [1, 1, 923]
)
MPATS = (['man', 'max'] if ID[-1]=='D' else ['man', 'max'])
if ID[-1] == 'C':
    COLS = ['#e01a4f77', '#a2d2ff00'] #'#a2d2ff']
GENS = 5000
###############################################################################
# File ID
###############################################################################
OUT_PTH = './sims_out/{}_{}'.format(ID, GENS)
# OUT_PTH = '/Users/sanchez.hmsc/Documents/WorkSims/MGSurvE_Validations/{}_{}'.format(ID, GENS)
srv.makeFolder(OUT_PTH)
###############################################################################
# Load Files
###############################################################################
mPat = MPATS[0]
(logs, traces) = ([], [])
for mPat in MPATS:
    logFiles = sorted(glob(path.join(OUT_PTH, (FPAT+'LOG.csv').format(mPat))))
    logs.append([pd.read_csv(f) for f in logFiles])
mins = [np.array([fc['min'] for fc in log]) for log in logs]
minFits = [m[-1] for m in mins[0]]
minVal = min(minFits) 
minIdx = minFits. index(minVal)
print('{} @ {}'.format(minVal, minIdx+1))
###############################################################################
# Plot GA Evolution
###############################################################################
(XRAN, YRAN) = ((0, 5000), (0, 150))
(fig, ax) = plt.subplots(figsize=(25, 3))
for (ix, trc) in enumerate(mins):
    ax.plot(trc.T/SCAL[ix], color=COLS[ix], lw=1.25)
ax.set_xlim(0, GENS)
ax.set_ylim(YRAN[0], YRAN[1])
ax.hlines(np.arange(YRAN[0], YRAN[1]+25, 25), XRAN[0], XRAN[1], color='#00000033', lw=1, zorder=-10)
ax.vlines(np.arange(XRAN[0], XRAN[1]+20, 500), YRAN[0], YRAN[1], color='#00000033', lw=1, zorder=-10)
ax.set_xticks([])
ax.set_yticks([])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
fig.savefig(
    path.join(OUT_PTH, (FPAT[:-7]+'-GA.png')), 
    facecolor=None, bbox_inches='tight', transparent=True,
    pad_inches=0, dpi=350
)
###############################################################################
# Load Landscape
###############################################################################
lndFiles = sorted(glob(path.join(OUT_PTH, (FPAT+'TRP.pkl').format('man'))))
lnd = srv.loadLandscape(
    OUT_PTH, lndFiles[0].split('/')[-1].split('.')[0], 
    fExt='pkl'
)
# Traps Kernels ---------------------------------------------------------------
(XRAN, YRAN) = ((0, 50), (0, 1))
(fig, ax) = plt.subplots(1, 1, figsize=(25, 3), sharey=False)
(fig, ax) = srv.plotTrapsKernels(
    fig, ax, lnd, distRange=(XRAN[0], XRAN[1]), aspect=.125
)
ax.set_xticks([])
ax.set_yticks([])
ax.hlines(np.arange(YRAN[0], YRAN[1]+.25, .25), XRAN[0], XRAN[1], color='#00000055', lw=1, zorder=-10)
ax.vlines(np.arange(XRAN[0], XRAN[1]+20, 12.5), YRAN[0], YRAN[1], color='#00000055', lw=1, zorder=-10)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
fig.savefig(
    path.join(OUT_PTH, '{}D-{}_KER.png'.format(ID, AP)), 
    facecolor=None, bbox_inches='tight', 
    pad_inches=0, dpi=300, transparent=True
)
plt.close('all')
###############################################################################
# Inspect Landscape
###############################################################################
PRINT_IDS = {'sites': False, 'traps': False}
(outer, itr, gen) = ('man', 1, 5000)
for outer in MPATS:
    logIx = MPATS.index(outer)
    if ID == 'YKND':
        itrsFitsPos = [aux.getBestTraps(l) for l in logs[logIx]]
        (fitVal, trpPos) = sorted(itrsFitsPos, key=lambda x: x[0])[0]
        fitsFun = aux.switchFunction(MPATS[logIx])
        trapsCoords = srv.chromosomeIDtoXY(trpPos, lnd.pointID, lnd.pointCoords).T
    else:
        itrsFitsPos = [aux.getBestTraps(l, discrete=False) for l in logs[logIx]]
        (fitVal, trpPos) = sorted(itrsFitsPos, key=lambda x: x[0])[0]
        fitsFun = aux.switchFunction(MPATS[logIx])
        trapsCoords = np.reshape(trpPos, (-1, 2)).T
    trapsLocs = pd.DataFrame(
        np.vstack([trapsCoords, lnd.trapsTypes, lnd.trapsFixed]).T, 
        columns=['lon', 'lat', 't', 'f']
    )
    trapsLocs['t'] = trapsLocs['t'].astype('int64')
    trapsLocs['f'] = trapsLocs['f'].astype('int64')
    # Update landscape --------------------------------------------------------
    lnd.updateTraps(trapsLocs, lnd.trapsKernels)
    if ID == 'YKND':
        fitness = srv.calcDiscreteFitness(
            trpPos, lnd, optimFunction=srv.getDaysTillTrapped,
            optimFunctionArgs={'inner': np.sum, 'outer': fitsFun}
        )[0]/SCAL[logIx]
    else:
        fitness = srv.calcFitness(
            trpPos, lnd, optimFunction=srv.getDaysTillTrapped,
            optimFunctionArgs={'inner': np.sum, 'outer': fitsFun}
        )[0]/SCAL[logIx]
    assert(np.isclose(fitVal, fitness*SCAL[logIx]))
    # Plot --------------------------------------------------------------------  
    (fig, ax) = (
        plt.figure(figsize=(15, 15)),
        plt.axes(projection=ccrs.PlateCarree())
    )
    lnd.plotSites(fig, ax, size=50)
    lnd.updateTrapsRadii([0.250, 0.125, 0.100])
    lnd.plotTraps(
        fig, ax, 
        zorders=(30, 25), transparencyHex='55', 
        proj=ccrs.PlateCarree()
    )
    srv.plotFitness(
        fig, ax, fitness, 
        fmt='{:.2f}', fontSize=25, color='#00000066', pos=(0.75, 0.10)
    )
    srv.plotClean(fig, ax, bbox=lnd.landLimits)
    if PRINT_IDS['sites']:
        for (ix, xy) in enumerate(lnd.pointCoords):
            ax.text(
                xy[0], xy[1], int(ix), 
                fontsize=3.5, zorder=10, ha='center', va='center_baseline'
            )
    if PRINT_IDS['traps']:
        for (ix, xy) in enumerate(lnd.trapsCoords):
            ax.text(
                xy[0], xy[1], int(ix), 
                fontsize=10, zorder=50, ha='center', va='center_baseline'
            )
    fig.savefig(
        path.join(OUT_PTH, (FPAT[:-1]+'.png').format(outer)), 
        facecolor='w', bbox_inches='tight', pad_inches=0.1, dpi=350,
        transparent=False
    )
    plt.close('all')
