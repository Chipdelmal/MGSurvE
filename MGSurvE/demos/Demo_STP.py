#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math
from vincenty import vincenty
import numpy as np
import pandas as pd
from os import path
from sys import argv
from copy import deepcopy
import matplotlib.pyplot as plt
from deap import base, creator, algorithms, tools
from compress_pickle import dump, load
from sklearn.preprocessing import normalize
import MGSurvE as srv
import warnings
warnings.filterwarnings('ignore', 'The iteration is not making good progress')


(ID, OUT_PTH) = (
    'STP', '/home/chipdelmal/Documents/WorkSims/MGSurvE_Benchmarks/STP/'
)
TRPS_NUM = 8
IX_SPLIT = 27
###############################################################################
# Load Pointset
###############################################################################
sites = pd.read_csv(path.join(OUT_PTH, 'stp_cluster_sites_pop_v5_fixed.csv'))
sites['t'] = [0] * sites.shape[0]
SAO_TOME_LL = sites.iloc[IX_SPLIT:]
SAO_bbox = (
    (min(SAO_TOME_LL['lon']), max(SAO_TOME_LL['lon'])),
    (min(SAO_TOME_LL['lat']), max(SAO_TOME_LL['lat']))
)
###############################################################################
# Load Migration Matrix
###############################################################################
migration = np.genfromtxt(
    path.join(OUT_PTH, 'kernel_cluster_v6a.csv'), delimiter=','
)
msplit = migration[IX_SPLIT:,IX_SPLIT:]
SAO_TOME_MIG = normalize(msplit, axis=1, norm='l2')
###############################################################################
# Defining Traps
###############################################################################
nullTraps = [0] * TRPS_NUM
(lonTrap, latTrap) = (
    np.random.uniform(SAO_bbox[0][0], SAO_bbox[0][1], TRPS_NUM),
    np.random.uniform(SAO_bbox[1][0], SAO_bbox[1][1], TRPS_NUM)
)
traps = pd.DataFrame({
    'lon': lonTrap, 'lat': latTrap,
    't': nullTraps, 'f': nullTraps
})
tKer = {0: {'kernel': srv.exponentialDecay, 'params': {'A': .5, 'b': 50}}}
###############################################################################
# Setting Landscape Up
###############################################################################
lnd = srv.Landscape(
    SAO_TOME_LL, migrationMatrix=SAO_TOME_MIG,
    traps=traps, trapsKernels=tKer,
    distanceFunction=vincenty
)
###############################################################################
# Plot Landscape
###############################################################################
(fig, ax) = plt.subplots(1, 1, figsize=(15, 15), sharey=False)
lnd.plotSites(fig, ax)
lnd.plotTraps(fig, ax)
srv.plotClean(fig, ax)