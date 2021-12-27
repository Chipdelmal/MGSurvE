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
IX_SPLIT = 27
###############################################################################
# Load Pointset
###############################################################################
sites = pd.read_csv(path.join(OUT_PTH, 'stp_cluster_sites_pop_v5_fixed.csv'))
sites['t'] = [0] * sites.shape[0]
SAO_TOME_LL = sites.iloc[IX_SPLIT:]
###############################################################################
# Load Migration Matrix
###############################################################################
migration = np.genfromtxt(
    path.join(OUT_PTH, 'kernel_cluster_v6a.csv'), 
    delimiter=','
)
msplit = migration[IX_SPLIT:,IX_SPLIT:]
SAO_TOME_MIG = normalize(msplit, axis=1, norm='l2')
###############################################################################
# Load Migration Matrix
###############################################################################
lnd = srv.Landscape(
    SAO_TOME_LL, migrationMatrix=SAO_TOME_MIG,
    distanceFunction=vincenty
)

(fig, ax) = plt.subplots(1, 1, figsize=(15, 15), sharey=False)
lnd.plotSites(fig, ax)
srv.plotClean(fig, ax)