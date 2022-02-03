#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math
from time import time
import numpy as np
import pandas as pd
from os import path
from sys import argv
from copy import deepcopy
import matplotlib.pyplot as plt
from deap import base, creator, algorithms, tools
from compress_pickle import dump, load
import MGSurvE as srv
import warnings
warnings.filterwarnings('ignore', 'The iteration is not making good progress')


(ID, OUT_PTH) = ('Grid', './sims_out/')
(ptsNum, ptsTypes) = (10, 3)
###############################################################################
# Defining Landscape
###############################################################################
# Mosquito movement kernel
mKer = {'params': [.075, 1.0e-10, math.inf], 'zeroInflation': .75}
# Pseudo-random landscape %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bbox = ((-80, 80), (-80, 80))
xy = srv.ptsRegularGrid(ptsNum, bbox).T
ptsNum = xy.shape[1]
# Generate landscape with one point-type ......................................
points_hom = pd.DataFrame({'x': xy[0], 'y': xy[1], 't': [0]*xy.shape[1]})
# Duplicate dataframe and replace to multiple point-types for heterogeneity ...
points_het = points_hom.copy()
points_het['t'] = np.random.choice(ptsTypes, ptsNum, p=[0.10, 0.70, 0.20])
msk = [
    [0.05, 0.70, 0.25],
    [0.30, 0.10, 0.60],
    [0.70, 0.10, 0.20],
]
###############################################################################
# Defining Traps
###############################################################################
nullTraps = [0, 0, 0]
traps = pd.DataFrame({
    'x': nullTraps, 'y': nullTraps, 't': [0, 0, 1], 'f': nullTraps
})
tKer = {
    0: {'kernel': srv.exponentialDecay, 'params': {'A': .5, 'b': .1}},
    1: {'kernel': srv.exponentialDecay, 'params': {'A': .25, 'b': .125}}
}
###############################################################################
# Setting Landscape Up
###############################################################################
lnd_hom = srv.Landscape(
    points_hom, 
    kernelParams=mKer, 
    traps=traps, trapsKernels=tKer
)
lnd_het = srv.Landscape(
    points_het, maskingMatrix=msk, 
    kernelParams=mKer,
    traps=traps, trapsKernels=tKer
)
###############################################################################
# Plot Landscapes
###############################################################################
# Homogeneous -----------------------------------------------------------------
bbox = lnd_hom.getBoundingBox()
trpMsk = srv.genFixedTrapsMask(lnd_hom.trapsFixed)
(fig, ax) = plt.subplots(1, 1, figsize=(15, 15), sharey=False)
lnd_hom.plotSites(fig, ax, size=200)
lnd_hom.plotMigrationNetwork(fig, ax, alphaMin=.6, lineWidth=25)
# lnd_hom.plotTraps(fig, ax)
srv.plotClean(fig, ax, frame=False)
fig.savefig(
    path.join(OUT_PTH, '{}_LND_HOM.png'.format(ID)), 
    facecolor='w', bbox_inches='tight', pad_inches=0, dpi=250
)
plt.close('all')
# Heterogeneous ---------------------------------------------------------------
bbox = lnd_het.getBoundingBox()
trpMsk = srv.genFixedTrapsMask(lnd_het.trapsFixed)
(fig, ax) = plt.subplots(1, 1, figsize=(15, 15), sharey=False)
lnd_het.plotSites(fig, ax, size=200)
lnd_het.plotMaskedMigrationNetwork(fig, ax, alphaMin=.6, lineWidth=25)
# lnd_het.plotTraps(fig, ax)
srv.plotClean(fig, ax, bbox=bbox, frame=False)
fig.savefig(
    path.join(OUT_PTH, '{}_LND_HET.png'.format(ID)), 
    facecolor='w', bbox_inches='tight', pad_inches=0, dpi=250
)
plt.close('all')
###############################################################################
# Dump Landscapes
###############################################################################
srv.dumpLandscape(lnd_hom, OUT_PTH, '{}_LND_HOM'.format(ID))
srv.dumpLandscape(lnd_het, OUT_PTH, '{}_LND_HET'.format(ID))