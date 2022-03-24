#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math
import numpy as np
import pandas as pd
import MGSurvE as srv
import matplotlib.pyplot as plt

(PT_OUT, DPI) = ('/home/chipdelmal/Documents/WorkSims/Mov/Demo/', 300)
PAD = 1
(minX, maxX) = (.5, 3.25)
(minY, maxY) = (0, 3)
bbox = ((-0.5, 3.5), (-1, 3.5))
###############################################################################
# XY Landscape with one point-type and one trap-type
###############################################################################
pts = pd.DataFrame({
    'x': [0, 3, 2.5, 1, 1.25, .5], 
    'y': [0, 2, 1, 0, .75, 1.25], 
    't': [0, 0, 0, 0, 0, 0],
})
msk = [
    [.2, .8],
    [.8, .2]
]
points = pd.DataFrame(pts, columns=('x', 'y', 't'))
# Traps info ------------------------------------------------------------------
traps = pd.DataFrame({
    'x': [1.75, 0.65], 
    'y': [0.25, 2.50], 
    't': [0, 1],
    'f': [0, 0]
})
tker = {
    0: {'kernel': srv.exponentialDecay, 'params': {'A': .30, 'b': 2}},
    1: {'kernel': srv.exponentialDecay, 'params': {'A': .50, 'b': 4}} 
}
# Land creation ---------------------------------------------------------------
lnd = srv.Landscape(points, maskingMatrix=msk, traps=traps, trapsKernels=tker)
###############################################################################
# Plot Traps
###############################################################################
kers = lnd.trapsKernels
srv.plotTrapsKernels(lnd)

