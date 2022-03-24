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
    0: {'kernel': srv.exponentialDecay, 'params': {'A': 1, 'b': 0.1}},
    1: {'kernel': srv.sigmoidDecay,     'params': {'A': 1.0, 'rate': 0.085, 'x0': 40}},
    2: {'kernel': srv.exponentialAttractiveness, 'params': {'A': 1, 'k': .025, 's': .2, 'gamma': .8, 'epsilon': 0}}, 
}
# Land creation ---------------------------------------------------------------
lnd = srv.Landscape(points, maskingMatrix=msk, traps=traps, trapsKernels=tker)
srv.plotTrapsKernels(lnd, distRange=(0, 100))
###############################################################################
# Plot Kernel
###############################################################################
lnd = srv.loadLandscape(
    '/Volumes/marshallShare/MGSurvE_Yorkeys/', 
    'YK2_06_TRP', fExt='pkl'
)
srv.plotTrapsKernels(lnd, distRange=(0, 100))
