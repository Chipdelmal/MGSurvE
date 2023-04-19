#!/usr/bin/env python
# -*- coding: utf-8 -*-

import gc
import time
import numpy as np
import pandas as pd
from os import path
from sys import argv
from glob import glob
import matplotlib
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import MGSurvE as srv
import auxiliary as aux
from PIL import Image
# matplotlib.use('agg')
import warnings
warnings.filterwarnings("ignore")


if srv.isNotebook():
    (ID, AP, RID) = ('YKND', 'man', '01')
else:
    (ID, AP, RID) = argv[1:]
(GENS, RID, FPAT) = (5000, int(RID), ID+'-{}_16*')
###############################################################################
# File ID
###############################################################################
FNAME = FPAT.format(AP)[:-1]+'-{:02d}_'.format(RID)
I_PTH = '/home/chipdelmal/Documents/WorkSims/MGSurvE_Validations/YKND_5000'
O_PTH = path.join(I_PTH, 'video')
srv.makeFolder(O_PTH)
###############################################################################
# Load Files
###############################################################################
(log, lnd) = (
    pd.read_csv(path.join(I_PTH, FNAME+'LOG.csv')),
    srv.loadLandscape(I_PTH, FNAME+'TRP', fExt='pkl')
)
###############################################################################
# Plot Clean Landscape
###############################################################################
(PROJ, FIGS, PAD, DPI) = (ccrs.PlateCarree(), (15, 15), 0, 350)
(fig, ax) = (plt.figure(figsize=FIGS), plt.axes(projection=PROJ))
lnd.plotSites(fig, ax, size=50)
srv.plotClean(fig, ax, bbox=lnd.landLimits)
fig.savefig(
    path.join(O_PTH, FNAME+'CLN.png'), 
    transparent=False, facecolor='w',
    bbox_inches='tight', pad_inches=PAD, dpi=DPI
)
###############################################################################
# Plot Optimization
###############################################################################
fitFun = aux.switchFunction(AP)

gen = 1000

trpEntry = log.iloc[gen]['traps']
if ID == 'YKND':
    trpPos = aux.idStringToArray(trpEntry, discrete=True)
    trpCds = srv.chromosomeIDtoXY(trpPos, lnd.pointID, lnd.pointCoords).T
# Assemble and update traps ---------------------------------------------------
trapsLocs = pd.DataFrame(
    np.vstack([trpCds, lnd.trapsTypes, lnd.trapsFixed]).T, 
    columns=['lon', 'lat', 't', 'f']
)
trapsLocs['t'] = trapsLocs['t'].astype('int64')
trapsLocs['f'] = trapsLocs['f'].astype('int64')
lnd.updateTraps(trapsLocs, lnd.trapsKernels)
# Calculate new fitness -------------------------------------------------------
if ID == 'YKND':
    fitness = srv.calcDiscreteFitness(
        trpPos, lnd, optimFunction=srv.getDaysTillTrapped,
        optimFunctionArgs={'inner': np.sum, 'outer': fitFun}
    )[0]
    
# Plot ------------------------------------------------------------------------  
(fig, ax) = (plt.figure(figsize=FIGS), plt.axes(projection=PROJ))
lnd.plotSites(fig, ax, size=50)
lnd.updateTrapsRadii([0.250, 0.125, 0.100])
lnd.plotTraps(
    fig, ax, 
    zorders=(30, 25), transparencyHex='55', 
    latlon=True, proj=PROJ
)
srv.plotFitness(
    fig, ax, fitness, 
    fmt='{:.5f}', fontSize=20, color='#00000066', pos=(0.75, 0.10)
)
srv.plotClean(fig, ax, bbox=lnd.landLimits)