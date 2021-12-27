#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math
from math import log
import numpy as np
import pandas as pd
from os import path
from sys import argv
from copy import deepcopy
import matplotlib.pyplot as plt
from deap import base, creator, algorithms, tools
from compress_pickle import dump, load
import MGSurvE as srv
from plots import plotDirectedNetwork
import MGSurvE.constants as cst
import warnings
import networkx as nx
from sklearn.preprocessing import normalize
import numpy as np
warnings.filterwarnings('ignore', 'The iteration is not making good progress')

###############################################################################
# Defining Landscape
###############################################################################
(ptsNum, PTS_TYPE) = (200, 3)
bbox = ((-225, 225), (-175, 175))
xy = srv.ptsRandUniform(ptsNum, bbox).T
mKer = {'params': [.075, 1.0e-10, math.inf], 'zeroInflation': .75}
pType = np.random.choice(PTS_TYPE, xy.shape[1])
points = pd.DataFrame({'x': xy[0], 'y': xy[1], 't': pType})
###############################################################################
# Defining Traps
###############################################################################
traps = pd.DataFrame({
    'x': [0, 10, 0, 100], 
    'y': [0, 20, -150, 5],
    't': [0, 1, 0, 1], 
    'f': [0, 0, 1, 1]
})
tKer = {
    0: {'kernel': srv.exponentialDecay, 'params': {'A': .5, 'b': .1}},
    1: {'kernel': srv.exponentialDecay, 'params': {'A': .25, 'b': .05}}
}
###############################################################################
# Setting Landscape Up
###############################################################################
lnd = srv.Landscape(
    points, kernelParams=mKer,
    traps=traps, trapsKernels=tKer
)
bbox = lnd.getBoundingBox()
trpMsk = srv.genFixedTrapsMask(lnd.trapsFixed)
(fig, ax) = plt.subplots(1, 1, figsize=(15, 15), sharey=False)
lnd.plotSites(fig, ax, size=100)
lnd.plotMigrationNetwork(fig, ax, alphaMin=.6, lineWidth=25)
lnd.plotTraps(fig, ax)
lnd.plotDirectedNetwork(fig, ax)
srv.plotClean(fig, ax, frame=False)
fig.savefig(
    path.join('NET_DEV.png'), 
    facecolor='w', bbox_inches='tight', pad_inches=0.1, dpi=300
)