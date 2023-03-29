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


if srv.isNotebook():
    (ID, AP, RID) = ('YKND', 'man', '01')
else:
    (ID, AP, RID) = argv[1:]
RID = int(RID)
(FPAT, MPATS) = ('YKND-{}_16*', ['sum', 'man', 'max'])
###############################################################################
# File ID
###############################################################################
GENS = 1000
OUT_PTH = './sims_out/'
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
COLS = ['#072ac8', '#f72585', '#8338ec']
SCAL = [923, 1, 1]
(fig, ax) = plt.subplots(figsize=(20, 8))
for (ix, trc) in enumerate(mins):
    ax.plot(trc.T/SCAL[ix], color=COLS[ix]+'88', lw=0.5)
ax.set_xlim(0, GENS)
ax.set_ylim(0, 150)
###############################################################################
# Load Landscape
###############################################################################
lndFiles = sorted(glob(path.join(OUT_PTH, (FPAT+'TRP.pkl').format('sum'))))
lnd = srv.loadLandscape(
    OUT_PTH, lndFiles[0].split('/')[-1].split('.')[0], 
    fExt='pkl'
)

probe = [500, 501, 527, 213, 688, 243, 531, 449, 703, 585, 115, 131, 212, 555, 101, 460]

gfend = 0
probe = np.array([int(i) for i in logs[0][0]['traps'].iloc[GENS-gfend][1:-1].split(',')])
trapsCoords = srv.chromosomeIDtoXY(probe, lnd.pointID, lnd.pointCoords).T
trapsLocs = pd.DataFrame(
    np.vstack([trapsCoords, lnd.trapsTypes, lnd.trapsFixed]).T, 
    columns=['lon', 'lat', 't', 'f']
)
trapsLocs['t']=trapsLocs['t'].astype('int64')
trapsLocs['f']=trapsLocs['f'].astype('int64')
lnd.updateTraps(trapsLocs, lnd.trapsKernels)

(fig, ax) = (
    plt.figure(figsize=(15, 15)),
    plt.axes(projection=ccrs.PlateCarree())
)
lnd.plotSites(fig, ax, size=50)
lnd.plotTraps(fig, ax, zorders=(30, 25))
# srv.plotFitness(fig, ax, min(logbook['min']), fmt='{:.5f}', fontSize=100)
srv.plotClean(fig, ax, bbox=lnd.landLimits)