#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math
import numpy as np
import pandas as pd
import MGSurvE as srv
import matplotlib.pyplot as plt

###############################################################################
# XY Landscape with one point-type and one trap-type
###############################################################################
pts = pd.DataFrame({'x': [0, 0], 'y': [0, 0], 't': [0, 0]})
msk = [[.2, .8], [.8, .2]]
points = pd.DataFrame(pts, columns=('x', 'y', 't'))
# Traps info ------------------------------------------------------------------
traps = pd.DataFrame({'x': [0, 0], 'y': [0, 0], 't': [0, 0], 'f': [0, 0]})
tKer = {
    0: {
        'kernel': srv.exponentialAttractiveness,
        'params': {'A': 1, 'k': .01, 's': .3, 'gamma': .975, 'epsilon': 0}
    },
    1: {
        'kernel': srv.exponentialDecay, 
        'params': {'A': 1, 'b': 0.045}
    },
    2: {
        'kernel': srv.sigmoidDecay,     
        'params': {'A': 1, 'rate': .175, 'x0': 25}
    }
}
TCOL = {
    0: '#f7258515', 1: '#fe5f5515', 2: '#5ddeb125', 
    3: '#f038ff15', 4: '#e2ef7015', 5: '#9381ff15', 
}
# Land creation ---------------------------------------------------------------
lnd = srv.Landscape(points, maskingMatrix=msk, traps=traps, trapsKernels=tKer)
(fig, ax) = plt.subplots(1, 1, figsize=(15, 15), sharey=False)
(fig, ax) = srv.plotTrapsKernels(
    fig, ax, lnd, 
    colors=TCOL, distRange=(0, 100), aspect=.1
)
fig.savefig(
    './kernels.png', 
    facecolor='w', bbox_inches='tight', pad_inches=0.1, dpi=300
)
###############################################################################
# Plot Kernel
###############################################################################
# lnd = srv.loadLandscape(
#     '/Volumes/marshallShare/MGSurvE_Yorkeys/', 
#     'YK2_06_TRP', fExt='pkl'
# )
# srv.plotTrapsKernels(lnd, distRange=(0, 100))
###############################################################################
# Kernel
###############################################################################
lnd = srv.Landscape(
    points, 
    maskingMatrix=msk, traps=traps, trapsKernels=tker,
    trapsRadii=[1, .1]
)
lnd.trapsKernels