#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math
import pandas as pd
import MGSurvE as srv
import matplotlib.pyplot as plt


OUT_PTH = './scratch/'
###############################################################################
# XY Landscape with one point-type and one trap-type
###############################################################################
pts = (
    (0.0, 0.0, 0), 
    (2.0, 0.5, 0), 
    (2.5, 1.5, 0),
)
points = pd.DataFrame(pts, columns=('x', 'y', 't'))
# Traps info ------------------------------------------------------------------
trp = (
    (2.5, 0.75, 0, 0),
    (0.0, 0.50, 0, 0)
)
traps = pd.DataFrame(trp, columns=('x', 'y', 't', 'f'))
tKernels = {0: {'kernel': srv.exponentialDecay, 'params': {'A': 0.5, 'b': 3}}}
# Land creation ---------------------------------------------------------------
lnd = srv.Landscape(points, traps=traps, trapsKernels=tKernels)
lnd.calcFundamentalMatrix()
lnd.getDaysTillTrapped()
# Plotting landscape ----------------------------------------------------------
(fig, ax) = plt.subplots(1, 2, figsize=(15, 15), sharey=False)
lnd.plotSites(fig, ax[0])
lnd.plotMigrationNetwork(fig, ax[0])
lnd.plotTraps(fig, ax[0])
lnd.plotTrapsNetwork(fig, ax[0])
srv.plotMatrix(fig, ax[1], lnd.trapsMigration, lnd.trapsNumber)
[srv.plotClean(fig, i, frame=False) for i in ax]
fig.savefig(
    './demo_basicLandscape.png', facecolor='w',
    bbox_inches='tight', pad_inches=0, dpi=300
)
###############################################################################
# Updating traps
###############################################################################
traps = pd.DataFrame({
    'x': [0.5, 3.0, 2.0], 
    'y': [0.0, 0.0, 2.0], 
    't': [0, 1, 0],
    'f': [1, 0, 0]
})
tKernels = {
    0: {'kernel': srv.exponentialDecay, 'params': {'A': .30, 'b': 2}},
    1: {'kernel': srv.exponentialDecay, 'params': {'A': .50, 'b': 1}} 
}
lnd.updateTraps(traps, tKernels)
# Plotting updated traps ------------------------------------------------------
(fig, ax) = plt.subplots(1, 2, figsize=(15, 15), sharey=False)
lnd.plotSites(fig, ax[0])
lnd.plotMigrationNetwork(fig, ax[0])
lnd.plotTraps(fig, ax[0])
lnd.plotTrapsNetwork(fig, ax[0])
srv.plotMatrix(fig, ax[1], lnd.trapsMigration, lnd.trapsNumber)
[srv.plotClean(fig, i, frame=False) for i in ax]
fig.savefig(
    path.join(OUT_PTH, 'demo_updatedLandscape.png'), facecolor='w',
    bbox_inches='tight', pad_inches=0, dpi=300
)
