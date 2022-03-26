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
tker = {
    0: {'kernel': srv.exponentialDecay, 'params': {'A': 1, 'b': 0.045}},
    1: {'kernel': srv.sigmoidDecay,     'params': {'A': 1, 'rate': .175, 'x0': 30}},
    2: {'kernel': srv.exponentialAttractiveness, 'params': {'A': 1, 'k': .01, 's': .3, 'gamma': .975, 'epsilon': 0}}, 
}
# Land creation ---------------------------------------------------------------
lnd = srv.Landscape(points, maskingMatrix=msk, traps=traps, trapsKernels=tker)
(fig, ax) = plt.subplots(1, 1, figsize=(15, 15), sharey=False)
(fig, ax) = srv.plotTrapsKernels(fig, ax, lnd, distRange=(0, 150))
# fig.savefig(
#     './kernels.png', 
#     facecolor='w', bbox_inches='tight', pad_inches=0.1, dpi=300
# )
###############################################################################
# Plot Kernel
###############################################################################
# lnd = srv.loadLandscape(
#     '/Volumes/marshallShare/MGSurvE_Yorkeys/', 
#     'YK2_06_TRP', fExt='pkl'
# )
# srv.plotTrapsKernels(lnd, distRange=(0, 100))
