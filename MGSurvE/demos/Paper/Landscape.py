#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import math
import warnings
import numpy as np
import pandas as pd
from os import path
import matplotlib.pyplot as plt
import Constants as cst
import MGSurvE as srv
warnings.filterwarnings('ignore', 'The iteration is not making good progress')

if srv.isNotebook():
    (ID, OUT_PTH, ZIK) = ("Grid", "./sims_out/", "ZI")
else:
    (ID, OUT_PTH, ZIK) = (sys.argv[1], cst.out_pth, sys.argv[2])
###############################################################################
# Synthetic Landscape Selector
###############################################################################
bbox = cst.bbox
if ID == 'Grid':
    (ptsNum, ptsTypes) = (int(math.sqrt(cst.ptsNum)), len(cst.pTypesProb))
    xy = srv.ptsRegularGrid(ptsNum, (bbox[0], bbox[0])).T
elif ID == 'Uniform':
    (ptsNum, ptsTypes) = (cst.ptsNum, len(cst.pTypesProb))
    xy = srv.ptsRandUniform(ptsNum, bbox).T
elif ID == 'Ring':
    (ptsNum, radii, ptsTypes) = (
        cst.ptsNum, (cst.bbox[1][0], int(cst.bbox[1][0]*0.6)), len(cst.pTypesProb)
    )
    xy = srv.ptsDonut(ptsNum, radii).T
elif ID == 'Poisson':
    (ptsNum, clsNum, radius, ptsTypes) = (
        cst.ptsNum, cst.clsNum, cst.clsRad, len(cst.pTypesProb)
    )
    xy = srv.ptsPossion(
        ptsNum, cst.clsNum, radius,
        bbox=bbox, polygon=None
    ).T
elif ID == 'Circle':
    (ptsNum, radius, ptsTypes) = (
        cst.ptsNum, cst.bbox[0][1], len(cst.pTypesProb)
    )
    xy = srv.ptsRegularCircle(ptsNum, radius, solStart=20).T
ptsNum = xy.shape[1]
# Generate landscape with one point-type ......................................
points_hom = pd.DataFrame({'x': xy[0], 'y': xy[1], 't': [0]*ptsNum})
# Duplicate dataframe and replace to multiple point-types for heterogeneity ...
points_het = points_hom.copy()
points_het['t'] = np.random.choice(ptsTypes, ptsNum, p=cst.pTypesProb)
###############################################################################
# Defining Traps
###############################################################################
nullT= cst.nullTraps
traps = pd.DataFrame({
    'sid': nullT,
    'x': nullT, 'y': nullT, 
    't': cst.typeTraps , 'f': nullT
})
tKer = cst.tKer
###############################################################################
# Setting Landscape Up
###############################################################################
(mKer, msk) = (
    (cst.mKerZ if ZIK=='ZI' else cst.mKerN), 
    cst.msk
)
lnd_hom = srv.Landscape(
    points_hom, 
    kernelParams=mKer, 
    traps=traps, trapsKernels=tKer,
    trapsRadii=[.5, .375, .25]
)
lnd_het = srv.Landscape(
    points_het, maskingMatrix=msk, 
    kernelParams=mKer,
    traps=traps, trapsKernels=tKer,
    trapsRadii=[.5, .375, .25]
)
###############################################################################
# Plot Landscapes
###############################################################################
ID = f'{ZIK}-{ID}'
# Homogeneous -----------------------------------------------------------------
bbox = lnd_hom.getBoundingBox()
trpMsk = srv.genFixedTrapsMask(lnd_hom.trapsFixed)
(fig, ax) = plt.subplots(1, 1, figsize=(15, 15), sharey=False)
lnd_hom.plotSites(fig, ax, size=200)
lnd_hom.plotMigrationNetwork(fig, ax, alphaMin=.5, lineWidth=50)
# lnd_hom.plotTraps(fig, ax)
srv.plotClean(fig, ax, frame=False, bbox=bbox, pad=cst.pad_i)
fig.savefig(
    path.join(OUT_PTH, '{}_LND_HOM.png'.format(ID)), 
    facecolor='w', bbox_inches='tight', pad_inches=cst.pad, dpi=cst.dpi
)
plt.close('all')
# Traps Kernels ---------------------------------------------------------------
(fig, ax) = plt.subplots(1, 1, figsize=(15, 5), sharey=False)
(fig, ax) = srv.plotTrapsKernels(fig, ax, lnd_hom, distRange=(0, 100), aspect=.175)
ax.set_xlabel("Distance (m)")
ax.set_ylabel("Attractiveness")
fig.savefig(
    path.join(OUT_PTH, '{}_{:02d}_KER.png'.format(ID, len(cst.typeTraps))), 
    facecolor='w', bbox_inches='tight', pad_inches=0.1, dpi=300
)
plt.close('all')
# Heterogeneous ---------------------------------------------------------------
bbox = lnd_het.getBoundingBox()
trpMsk = srv.genFixedTrapsMask(lnd_het.trapsFixed)
(fig, ax) = plt.subplots(1, 1, figsize=(15, 15), sharey=False)
lnd_het.plotSites(fig, ax, size=200)
lnd_het.plotMaskedMigrationNetwork(fig, ax, alphaMin=.5, lineWidth=50)
# lnd_het.plotTraps(fig, ax)
srv.plotClean(fig, ax, frame=False, bbox=bbox, pad=cst.pad_i)
fig.savefig(
    path.join(OUT_PTH, '{}_LND_HET.png'.format(ID)), 
    facecolor='w', bbox_inches='tight', pad_inches=cst.pad, dpi=cst.dpi
)
plt.close('all')
###############################################################################
# Dump Landscapes
###############################################################################
srv.dumpLandscape(lnd_hom, OUT_PTH, '{}_LND_HOM'.format(ID))
srv.dumpLandscape(lnd_het, OUT_PTH, '{}_LND_HET'.format(ID))
