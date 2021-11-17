#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math
import numpy as np
import pandas as pd
import MGSurvE as srv
import matplotlib.pyplot as plt
from sklearn.preprocessing import normalize

###############################################################################
# XY Landscape
###############################################################################
pts = [
    [0.00, 0.00, 0], 
    [0.25, 0.50, 1], 
    [2.5, 0.15, 0],
    [5, 0.1, 0],
    [3, 3, 0]
]
points = pd.DataFrame(pts, columns=['x', 'y', 't'])
msk = [
    [.6, .4],
    [.3, .7]
]
# Traps info ------------------------------------------------------------------
trp = [
    [5, 1, 0],
    [0, .5, 1]
]
traps = pd.DataFrame(trp, columns=['x', 'y', 't'])
tker = {
    0: {'kernel': srv.exponentialDecay, 'params': {'A': 0.5, 'b': 2}},
    1: {'kernel': srv.exponentialDecay, 'params': {'A': 0.5, 'b': 2}} 
}
# Land tests ------------------------------------------------------------------
lnd = srv.Landscape(
    points, maskingMatrix=msk, traps=traps, trapsKernels=tker
)
lnd.distanceMatrix
lnd.migrationMatrix
lnd.maskedMigration
lnd.trapsMigration


(fig, ax) = plt.subplots(figsize=(15, 15))
srv.plotMatrix(fig, ax, lnd.trapsMigration, lnd.trapsNumber)
srv.plotClean(fig, ax)


(fig, ax) = plt.subplots(figsize=(15, 15))
lnd.plotSites(fig, ax)
lnd.plotMigrationNetwork(fig, ax)
lnd.plotTraps(fig, ax)
lnd.plotTrapsNetwork(fig, ax)
srv.plotClean(fig, ax)



trapsRadii=[.5, .01]
lnd.updateTrapsRadii(trapsRadii)
###############################################################################
# Active dev
###############################################################################
trapsCoords = lnd.trapsCoords
trapsDistances = lnd.trapsDistances
trapsTypes = lnd.trapsTypes
trapsKernels = lnd.trapsKernels
pointCoords = lnd.pointCoords
trapsNumber = lnd.trapsNumber
pointNumber = lnd.pointNumber

traps = pd.DataFrame({
    'x': [0, 3],
    'y': [1, 0],
    't': [0, 2]
})
tker = {
    0: {'kernel': srv.exponentialDecay, 'params': {'A': .25, 'b': 1.5}},
    2: {'kernel': srv.exponentialDecay, 'params': {'A': .25, 'b': 2}} 
}
lnd.updateTraps(traps, tker)

trapsTypes = lnd.trapsTypes
trapsCoords = lnd.trapsCoords
# Plots tests -----------------------------------------------------------------
(fig, ax) = plt.subplots(figsize=(15, 15))
(fig, ax) = srv.plotSites(
    fig, ax, 
    lnd.pointCoords, lnd.pointTypes,
    size=500, edgecolors='w', linewidths=1.25,
    zorder=5
)
(fig, ax) = srv.plotMigrationNetwork(
    fig, ax, 
    lnd.maskedMigration, lnd.pointCoords, lnd.pointCoords,
    lineWidth=20, alphaMin=.5, alphaAmplitude=2.5, 
    zorder=0
)
(fig, ax) = srv.plotTraps(
    fig, ax, 
    trapsCoords, trapsTypes, lnd.trapsKernels,
    colors=srv.TRP_COLS
)
srv.plotTrapsNetwork(
    fig, ax,
    lnd.trapsMigration, lnd.trapsCoords, lnd.pointCoords
)
ax.set_aspect('equal')



(fig, ax) = plt.subplots(figsize=(15, 15))
srv.plotMatrix(fig, ax, lnd.trapsMigration)


###############################################################################
# Geo Landscape
###############################################################################
# stp = [
#     [7.353119999999999656e+00,1.598879999999999857e+00],
#     [7.377180000000000071e+00,1.620500000000000052e+00],
#     [7.379510000000000680e+00,1.678670000000000107e+00],
#     [7.380225000000000257e+00,1.652884999999999938e+00],
#     [7.384980000000000544e+00,1.670090000000000074e+00],
#     [7.388799999999999812e+00,1.632660000000000000e+00],
#     [7.398530000000000051e+00,1.635620000000000074e+00],
#     [7.399319999999999453e+00,1.546259999999999968e+00],
#     [7.410453332999999532e+00,1.671453333000000097e+00],
#     [7.420322221000000162e+00,1.645302658000000084e+00],
#     [7.413835000000000619e+00,1.620854999999999935e+00],
#     [7.417660000000000586e+00,1.607199999999999962e+00],
#     [7.431000000000000050e+00,1.624659999999999993e+00],
#     [7.437073332999999842e+00,1.681186666999999968e+00],
#     [7.441730000000000622e+00,1.668490000000000251e+00],
#     [7.453380000000000116e+00,1.681509999999999838e+00],
#     [7.457160000000000899e+00,1.624570000000000070e+00],
#     [6.465444999999999887e+00,2.474349999999999883e-01],
#     [6.476387653999999827e+00,2.628326810000000124e-01],
#     [6.475876315999999910e+00,2.176569969999999909e-01],
#     [6.476329999999999920e+00,2.349699999999999844e-01],
#     [6.486159999999999926e+00,2.854099999999999970e-01],
#     [6.491584999999999717e+00,3.015849999999999920e-01],
#     [6.494855000000000267e+00,1.272150000000000225e-01],
#     [6.509799425999999833e+00,2.971222589999999997e-01],
#     [6.499249999999999972e+00,2.229699999999999738e-01],
#     [6.500459999999999461e+00,3.151400000000000312e-01],
#     [6.504662913999999851e+00,2.782042809999999977e-01],
#     [6.506962601000000568e+00,2.650615749999999937e-01],
#     [6.517344999999999722e+00,3.294324999999999615e-01],
#     [6.514495608999999909e+00,2.184257139999999930e-01],
#     [6.515209999999999724e+00,1.022899999999999920e-01]
# ]
# points = pd.DataFrame(stp, columns=['lon', 'lat'])
# lnd = srv.Landscape(points)
# lnd.distanceMatrix
# lnd.migrationMatrix