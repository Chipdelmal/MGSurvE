#!/usr/bin/env python
# -*- coding: utf-8 -*-

from glob import glob
import numpy as np
import pandas as pd
from os import path
from sys import argv
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import MGSurvE as srv
import auxiliary as aux
from copy import deepcopy
plt.rcParams['axes.facecolor']='white'
plt.rcParams['savefig.facecolor']='white'


if srv.isNotebook():
    (ID, AP, RID) = ('YKND', 'man', '01')
else:
    (ID, AP, RID) = argv[1:]
RID = int(RID)
FPAT = ID+'-{}_16*'
(MPATS, COLS, SCAL) =  (
    ['man', ], # 'sum' 'max'],
    ['#072ac8', '#f72585', '#8338ec'],
    [1, 923, 1]
)
GENS = 5000
###############################################################################
# File ID
###############################################################################
OUT_PTH = '/home/chipdelmal/Documents/WorkSims/MGSurvE_Validations/YKND_5000'
# OUT_PTH = '/Users/sanchez.hmsc/Documents/GitHub/MGSurvE/MGSurvE/demos/YKN_Validations/sims_out'
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
###############################################################################
# Plot GA Evolution
###############################################################################
(fig, ax) = plt.subplots(figsize=(20, 6))
for (ix, trc) in enumerate(mins):
    ax.plot(trc.T/SCAL[ix], color=COLS[ix]+'66', lw=2)
ax.set_xlim(0, GENS)
ax.set_ylim(50, 100)
fig.savefig(
    path.join(OUT_PTH, (FPAT[:-7]+'-GA.png')), 
    facecolor='w', bbox_inches='tight', pad_inches=0.1, dpi=300
)
###############################################################################
# Load Landscape
###############################################################################
lndFiles = sorted(glob(path.join(OUT_PTH, (FPAT+'TRP.pkl').format('man'))))
lnd = srv.loadLandscape(
    OUT_PTH, lndFiles[0].split('/')[-1].split('.')[0], 
    fExt='pkl'
)
###############################################################################
# Inspect Landscape
#
# probe = [500, 501, 527, 213, 688, 243, 531, 449, 703, 585, 115, 131, 212, 555, 101, 460]
###############################################################################
PRINT_IDS = {'sites': False, 'traps': True}
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
    # lnd.plotMigrationNetwork(
    #     fig, ax, 
    #     lineWidth=20, alphaMin=.2, alphaAmplitude=15
    # )
    lnd.plotTraps(fig, ax, zorders=(30, 25), transparencyHex='55')
    srv.plotFitness(
        fig, ax, fitness, 
        fmt='{:.5f}', fontSize=20, color='#00000066', pos=(0.75, 0.10)
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
        facecolor='w', bbox_inches='tight', pad_inches=0.1, dpi=300,
        transparent=False
    )
###############################################################################
# Fine-grained optimization for discrete case
###############################################################################
# BRUTE_FORCE_TRAP_ID = 11
# (outer, itr, gen) = ('man', 1, 5000)
# lndBst = deepcopy(lnd)
# if ID == 'YKND':
#     itrsFitsPos = [aux.getBestTraps(l) for l in logs[logIx]]
#     (fitVal, trpPos) = sorted(itrsFitsPos, key=lambda x: x[0])[0]
#     bGA = {'fitVal': fitVal, 'trpPos': np.copy(trpPos)}
#     fitsFun = aux.switchFunction(MPATS[logIx])
#     for i in range(lnd.pointCoords.shape[0]):
#         trpPos[BRUTE_FORCE_TRAP_ID] = i
#         trapsCoords = srv.chromosomeIDtoXY(trpPos, lnd.pointID, lnd.pointCoords).T
#         trapsLocs = pd.DataFrame(
#             np.vstack([trapsCoords, lnd.trapsTypes, lnd.trapsFixed]).T, 
#             columns=['lon', 'lat', 't', 'f']
#         )
#         trapsLocs['t'] = trapsLocs['t'].astype('int64')
#         trapsLocs['f'] = trapsLocs['f'].astype('int64')
#         lnd.updateTraps(trapsLocs, lnd.trapsKernels)
#         fitness = srv.calcDiscreteFitness(
#             trpPos, lnd, optimFunction=srv.getDaysTillTrapped,
#             optimFunctionArgs={'inner': np.sum, 'outer': fitsFun}
#         )[0]/SCAL[logIx]
#         if fitness <= bGA['fitVal']:
#             print("Better position found! (fitness: {})".format(fitness))
#             bGA['fitVal'] = fitness
#             lndBst = deepcopy(lnd)
#     # Plot --------------------------------------------------------------------
#     (fig, ax) = (
#         plt.figure(figsize=(15, 15)),
#         plt.axes(projection=ccrs.PlateCarree())
#     )
#     lndBst.plotSites(fig, ax, size=50)
#     lndBst.plotTraps(fig, ax, zorders=(30, 25), transparencyHex='55')
#     srv.plotFitness(
#         fig, ax, bGA['fitVal'], 
#         fmt='{:.5f}', fontSize=20, color='#00000066', pos=(0.75, 0.10)
#     )
#     srv.plotClean(fig, ax, bbox=lndBst.landLimits)
#     fig.savefig(
#         path.join(OUT_PTH, (FPAT[:-1]+'-BRUTE.png').format(outer)), 
#         facecolor='w', bbox_inches='tight', pad_inches=0.1, dpi=300,
#         transparent=False
#     )