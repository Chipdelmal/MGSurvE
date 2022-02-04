#!/usr/bin/env python
# -*- coding: utf-8 -*-

import warnings
import numpy as np
import pandas as pd
from os import path
import matplotlib.pyplot as plt
import Constants as cst
import MGSurvE as srv
warnings.filterwarnings('ignore', 'The iteration is not making good progress')


(ID, OUT_PTH) = ('Poisson', './sims_out/')
(ptsNum, clsNum, radius, ptsTypes) = (
    cst.ptsNum, cst.clsNum, cst.clsRad, len(cst.pTypesProb)
)
###############################################################################
# Defining Landscape
###############################################################################
# Mosquito movement kernel
mKer = cst.mKer
# Pseudo-random landscape %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bbox = cst.bbox
xy = srv.clusterPossion(
    ptsNum, cst.clsNum, radius,
    bbox=bbox, polygon=None
).T
ptsNum = xy.shape[1]
# Generate landscape with one point-type ......................................
points_hom = pd.DataFrame({'x': xy[0], 'y': xy[1], 't': [0]*xy.shape[1]})
# Duplicate dataframe and replace to multiple point-types for heterogeneity ...
points_het = points_hom.copy()
points_het['t'] = np.random.choice(ptsTypes, ptsNum, p=cst.pTypesProb)
msk = cst.msk
###############################################################################
# Defining Traps
###############################################################################
nullTraps = cst.nullTraps
traps = pd.DataFrame({
    'x': nullTraps, 'y': nullTraps, 't': cst.typeTraps , 'f': nullTraps
})
tKer = cst.tKer
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
srv.plotClean(fig, ax, frame=False)
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