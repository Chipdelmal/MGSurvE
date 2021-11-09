#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import MGSurvE as srv
import matplotlib.pyplot as plt

pts = [
    [0.00, 0.00, 0], 
    [0.25, 0.50, 1], 
    [1.00, 0.15, 0]
]
points = pd.DataFrame(pts, columns=['x', 'y', 't'])

msk = [
    [.6, .4],
    [.3, .7]
]

lnd = srv.Landscape(points, maskingMatrix=msk)

lnd.distanceMatrix
lnd.migrationMatrix
lnd.maskedMigration

(fig, ax) = plt.subplots(figsize=(15, 15))
srv.plotSites(
    fig, ax, 
    lnd.pointCoords, lnd.pointTypes,
    size=500, edgecolors='w', linewidths=1.25,
    zorder=5
)
srv.plotNetwork(
    fig, ax, 
    lnd.maskedMigration, lnd.pointCoords,
    lineWidth=20, alphaMin=.5, alphaAmplitude=2.5, 
    zorder=0
)

