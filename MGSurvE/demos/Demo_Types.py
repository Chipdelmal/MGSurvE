#!/usr/bin/env python
# -*- coding: utf-8 -*-

from os import path
import pandas as pd
import MGSurvE as srv
import matplotlib.pyplot as plt


OUT_PTH = './scratch/'
###############################################################################
# XY Landscape with two point types
###############################################################################
pts = [
    [-4.0, 4.00, 0], 
    [0.25, 8.00, 1], 
    [5.00, 0.15, 0],
    [-1.0, 1.00, 0],
    [3.00, 3.00, 1]
]
points = pd.DataFrame(pts, columns=['x', 'y', 't'])
msk = [
    [0.05, 0.95],
    [0.95, 0.05]
]
###############################################################################
# XY Landscape with two trap types
###############################################################################
trp = [
    [5.00, 2.00, 1, 0],
    [-2.0, 2.00, 0, 0],
    [10.0, 0.00, 0, 1],
]
traps = pd.DataFrame(trp, columns=['x', 'y', 't', 'f'])
tker = {
    0: {'kernel': srv.exponentialDecay, 'params': {'A': 0.4, 'b': .5}},
    1: {'kernel': srv.sigmoidDecay,     'params': {'A': .6, 'rate': 1, 'x0': 0}}
}
###############################################################################
# Landscape Object
###############################################################################
lnd = srv.Landscape(
    points, maskingMatrix=msk, traps=traps, trapsKernels=tker
)
###############################################################################
# Plot
###############################################################################
(fig, ax) = plt.subplots(1, 2, figsize=(15, 15), sharey=False)
lnd.plotSites(fig, ax[0])
lnd.plotMaskedMigrationNetwork(fig, ax[0])
lnd.plotTraps(fig, ax[0])
lnd.plotTrapsNetwork(fig, ax[0])
srv.plotMatrix(fig, ax[1], lnd.trapsMigration, lnd.trapsNumber)
[srv.plotClean(fig, i, frame=False) for i in ax]
fig.savefig(
    path.join(OUT_PTH, 'demo_pointTypes.png'), facecolor='w',
    bbox_inches='tight', pad_inches=0, dpi=300
)