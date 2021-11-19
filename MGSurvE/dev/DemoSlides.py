#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import MGSurvE as srv
import matplotlib.pyplot as plt

(PT_OUT, DPI) = ('/home/chipdelmal/Documents/WorkSims/Mov/Demo/', 300)
PAD = 1
(minX, maxX) = (0, 3.25)
(minY, maxY) = (0, 3)
###############################################################################
# XY Landscape with one point-type and one trap-type
###############################################################################
pts = pd.DataFrame({
    'x': [0, 3, 2.5, 1, 1], 
    'y': [0, 2, 1, 0, 2.25], 
    't': [0, 0, 1, 0, 0],
})
points = pd.DataFrame(pts, columns=('x', 'y', 't'))
msk = [
    [.2, .8],
    [.8, .2]
]
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
# Fitness ---------------------------------------------------------------------
(tau, sitesN, trapsN) = (lnd.trapsMigration, lnd.pointNumber, lnd.trapsNumber)
fitFuns = {'outer': np.mean, 'inner': np.mean}
F = srv.getFundamentalMatrix(tau, sitesN, trapsN)
daysTillTrapped = np.apply_along_axis(fitFuns['inner'], 1, F)
fitness = fitFuns['outer'](daysTillTrapped)
###############################################################################
# Plotting landscape
###############################################################################
(fig, ax) = plt.subplots(1, 2, figsize=(15, 15), sharey=False)
lnd.plotSites(fig, ax[0])
lnd.plotMigrationNetwork(fig, ax[0])
srv.plotMatrix(fig, ax[1], lnd.migrationMatrix, vmax=.5)
[srv.plotClean(fig, i, frame=False) for i in ax]
fig.savefig(
    PT_OUT+'01.png', dpi=300, bbox_inches='tight', 
    pad_inches=0, transparent=False
)
# Plotting landscape ----------------------------------------------------------
(fig, ax) = plt.subplots(1, 2, figsize=(15, 15), sharey=False)
lnd.plotSites(fig, ax[0])
lnd.plotMaskedMigrationNetwork(fig, ax[0])
srv.plotMatrix(fig, ax[1], lnd.maskedMigration, vmax=.5)
[srv.plotClean(fig, i, frame=False) for i in ax]
fig.savefig(
    PT_OUT+'02.png', dpi=300, bbox_inches='tight', 
    pad_inches=0, transparent=False
)
# Plotting traps --------------------------------------------------------------
(fig, ax) = plt.subplots(1, 2, figsize=(15, 15), sharey=False)
lnd.plotSites(fig, ax[0])
lnd.plotMigrationNetwork(fig, ax[0])
lnd.plotTraps(fig, ax[0])
lnd.plotTrapsNetwork(fig, ax[0])
srv.plotMatrix(fig, ax[1], lnd.trapsMigration, lnd.trapsNumber)
[srv.plotClean(fig, i, frame=False) for i in ax]
fig.savefig(
    PT_OUT+'03.png', dpi=300, bbox_inches='tight', 
    pad_inches=0, transparent=False
)
# Plotting traps --------------------------------------------------------------
(fig, ax) = plt.subplots(1, 2, figsize=(15, 15), sharey=False)
lnd.plotSites(fig, ax[0])
lnd.plotMaskedMigrationNetwork(fig, ax[0])
lnd.plotTraps(fig, ax[0])
lnd.plotTrapsNetwork(fig, ax[0])
srv.plotMatrix(fig, ax[1], lnd.trapsMigration, lnd.trapsNumber)
[srv.plotClean(fig, i, frame=False) for i in ax]
fig.savefig(
    PT_OUT+'04.png', dpi=300, bbox_inches='tight', 
    pad_inches=0, transparent=False
)
###############################################################################
# Plotting traps locations
###############################################################################
(fig, ax) = plt.subplots(1, 1, figsize=(15, 15), sharey=False)
lnd.plotSites(fig, ax)
lnd.plotMaskedMigrationNetwork(fig, ax)
lnd.plotTraps(fig, ax)
lnd.plotTrapsNetwork(fig, ax)
srv.plotClean(fig, ax, frame=False)
ax.set_xlim(minX-PAD, maxX+PAD)
ax.set_ylim(minY-PAD, maxY+PAD)
ax.text(
    0.5, 0.5, '{:.2f}'.format(fitness),
    horizontalalignment='center', verticalalignment='center',
    fontsize=100, color='#00000011',
    transform=ax.transAxes, zorder=5
)
fig.savefig(
    PT_OUT+'05.png', dpi=300, bbox_inches='tight', 
    pad_inches=0, transparent=False
)
# New layout ------------------------------------------------------------------
traps = pd.DataFrame({
    'x': [3, 0.65], 
    'y': [0.25, 3], 
    't': [0, 1], 'f': [0, 0]
})
lnd.updateTraps(traps, trapsKernels=tker)
(tau, sitesN, trapsN) = (lnd.trapsMigration, lnd.pointNumber, lnd.trapsNumber)
F = srv.getFundamentalMatrix(tau, sitesN, trapsN)
daysTillTrapped = np.apply_along_axis(fitFuns['inner'], 1, F)
fitness = fitFuns['outer'](daysTillTrapped)
# Plot 
(fig, ax) = plt.subplots(1, 1, figsize=(15, 15), sharey=False)
lnd.plotSites(fig, ax)
lnd.plotMaskedMigrationNetwork(fig, ax)
lnd.plotTraps(fig, ax)
lnd.plotTrapsNetwork(fig, ax)
srv.plotClean(fig, ax, frame=False)    
ax.set_xlim(minX-PAD, maxX+PAD)
ax.set_ylim(minY-PAD, maxY+PAD)
ax.text(
    0.5, 0.5, '{:.2f}'.format(fitness),
    horizontalalignment='center', verticalalignment='center',
    fontsize=100, color='#00000011',
    transform=ax.transAxes, zorder=5
)
fig.savefig(
    PT_OUT+'06.png', dpi=300, bbox_inches='tight', 
    pad_inches=0, transparent=False
)
# New layout ------------------------------------------------------------------
traps = pd.DataFrame({
    'x': [2.5, .5], 
    'y': [1.5, 0], 
    't': [0, 1], 'f': [0, 0]
})
lnd.updateTraps(traps, trapsKernels=tker)
(tau, sitesN, trapsN) = (lnd.trapsMigration, lnd.pointNumber, lnd.trapsNumber)
F = srv.getFundamentalMatrix(tau, sitesN, trapsN)
daysTillTrapped = np.apply_along_axis(fitFuns['inner'], 1, F)
fitness = fitFuns['outer'](daysTillTrapped)
# Plot 
(fig, ax) = plt.subplots(1, 1, figsize=(15, 15), sharey=False)
lnd.plotSites(fig, ax)
lnd.plotMaskedMigrationNetwork(fig, ax)
lnd.plotTraps(fig, ax)
lnd.plotTrapsNetwork(fig, ax)
srv.plotClean(fig, ax, frame=False)
ax.set_xlim(minX-PAD, maxX+PAD)
ax.set_ylim(minY-PAD, maxY+PAD)
ax.text(
    0.5, 0.5, '{:.2f}'.format(fitness),
    horizontalalignment='center', verticalalignment='center',
    fontsize=100, color='#00000011',
    transform=ax.transAxes, zorder=5
)
fig.savefig(
    PT_OUT+'07.png', dpi=300, bbox_inches='tight', 
    pad_inches=0, transparent=False
)