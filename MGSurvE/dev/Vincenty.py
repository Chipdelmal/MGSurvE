#!/usr/bin/env python
# -*- coding: utf-8 -*-


import dill
import math
import numpy as np
import pandas as pd
from os import path
from sys import argv
from copy import deepcopy
import cartopy.crs as crs
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from haversine import haversine
from vincenty import vincenty
from compress_pickle import dump, load
from deap import base, creator, algorithms, tools
import MGSurvE as srv
import warnings
warnings.filterwarnings('ignore', 'The iteration is not making good progress')


(LND_PTH, OUT_PTH, ID, EXP) = (
    '/home/chipdelmal/Documents/WorkSims/MGSurvE_Yorkeys/LandOriginal/Yorkeys02.csv',
    '/home/chipdelmal/Documents/WorkSims/MGSurvE_Yorkeys/', 
    'YK2', '001'
)
GENS = 1000
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
# YK_CNTR = [i[0]+(i[1]-i[0])/2 for i in YK_LL]
# Movement Kernel -------------------------------------------------------------
mKer = {
    'kernelFunction': srv.zeroInflatedExponentialKernel,
    'kernelParams': {'params': srv.AEDES_EXP_PARAMS, 'zeroInflation': 0.0}
}
###############################################################################
# Defining Traps
###############################################################################
TRPS_NUM = 4
nullTraps = [0]*TRPS_NUM
traps = pd.DataFrame({
    'lon': [np.mean(YK_LL['lon'])]*TRPS_NUM, 'lat': [np.mean(YK_LL['lat'])]*TRPS_NUM,
    't': [0, 0, 0, 0], 'f': nullTraps
})
tKer = {
    0: {'kernel': srv.exponentialDecay, 'params': {'A': 1.0, 'b': 0.1               }},
    1: {'kernel': srv.sigmoidDecay,     'params': {'A': 1.0, 'rate': 0.2, 'x0': 25  }},
}
###############################################################################
# Setting Landscape Up
###############################################################################
lnd = srv.Landscape(
    YK_LL, 
    # distanceFunction=(lambda lat, lon: haversine(lat, lon, unit='m')),
    kernelFunction=mKer['kernelFunction'], kernelParams=mKer['kernelParams'],
    traps=traps, trapsKernels=tKer, trapsRadii=[.20, .25, .5],
    landLimits=YK_BBOX
)
bbox = lnd.getBoundingBox()
trpMsk = srv.genFixedTrapsMask(lnd.trapsFixed)
###############################################################################
# Plot Landscape
###############################################################################
(fig, ax) = (plt.figure(figsize=(15, 15)), plt.axes(projection=crs.PlateCarree()))
lnd.plotSites(fig, ax, size=75)
lnd.plotMigrationNetwork(
    fig, ax, 
    lineWidth=10, alphaMin=.1, alphaAmplitude=2.5,
)
lnd.plotTraps(fig, ax, zorders=(30, 25))
srv.plotClean(fig, ax, bbox=YK_BBOX)
fig.savefig(
    path.join(OUT_PTH, '{}{}_CLN.png'.format(OUT_PTH, ID, TRPS_NUM)), 
    facecolor='w', bbox_inches='tight', pad_inches=0.1, dpi=300
)
plt.close('all')
###############################################################################
# Debugging
###############################################################################
srv.dumpLandscape(lnd, OUT_PTH, '{}_{:02d}_TRP'.format(ID, TRPS_NUM), fExt='.pkl')




pth = path.join(OUT_PTH, '{}_{:02d}_TRP'.format(ID, TRPS_NUM))
with open(pth+'.pkl', 'wb') as f:
    pickle.dump(lnd, f)

with open(pth+'.pkl', 'rb') as f:
    loaded_obj = pickle.load(f)