#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math
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

(ID, OUT_PTH) = ('Ring', './sims_out/')
(ptsNum, radii, ptsTypes) = (100, (75, 100), 3)
###############################################################################
# Defining Landscape
###############################################################################
# Mosquito movement kernel
mKer = {'params': [.075, 1.0e-10, math.inf], 'zeroInflation': .75}
# Pseudo-random landscape %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xy = srv.ptsDonut(ptsNum, radii).T
# Generate landscape with one point-type ......................................
points_hom = pd.DataFrame({'x': xy[0], 'y': xy[1], 't': [0]*xy.shape[1]})
# Duplicate dataframe and replace to multiple point-types for heterogeneity ...
points_het = points_null.copy()
points_het['t'] = np.random.choice(ptsTypes, xy.shape[1])
###############################################################################
# Defining Traps
###############################################################################
nullTraps = [0, 0, 0, 0]
traps = pd.DataFrame({
    'x': nullTraps, 'y': nullTraps, 't': nullTraps, 'f': nullTraps
})
tKer = {0: {'kernel': srv.exponentialDecay, 'params': {'A': .5, 'b': .1}}}
###############################################################################
# Setting Landscape Up
###############################################################################
lnd_hom = srv.Landscape(
    points_hom, kernelParams=mKer, traps=traps, trapsKernels=tKer
)
lnd_het = srv.Landscape(
    points_het, kernelParams=mKer, traps=traps, trapsKernels=tKer
)
###############################################################################
# Plot Landscapes
###############################################################################
# bbox = lnd.getBoundingBox()
# trpMsk = srv.genFixedTrapsMask(lnd.trapsFixed)
# (fig, ax) = plt.subplots(1, 1, figsize=(15, 15), sharey=False)
# lnd.plotSites(fig, ax, size=100)
# lnd.plotMigrationNetwork(fig, ax, alphaMin=.6, lineWidth=25)
# lnd.plotTraps(fig, ax)
# srv.plotClean(fig, ax, frame=False)
###############################################################################
# Dump Landscapes
###############################################################################
srv.dumpLandscape(lnd_hom, OUT_PTH, '{}_LND_HOM'.format(ID))
srv.dumpLandscape(lnd_het, OUT_PTH, '{}_LND_HET'.format(ID))