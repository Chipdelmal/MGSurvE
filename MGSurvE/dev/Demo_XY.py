#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
import MGSurvE as srv
import matplotlib.pyplot as plt

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
    (2.5, 0.75, 0),
    (0, 0.5, 0)
)
traps = pd.DataFrame(trp, columns=('x', 'y', 't'))
tKernels = {0: {'kernel': srv.exponentialDecay, 'params': {'A': 0.5, 'b': 3}}}
# Land tests ------------------------------------------------------------------
lnd = srv.Landscape(points, traps=traps, trapsKernels=tKernels)
###############################################################################
# Plotting
###############################################################################
(fig, ax) = plt.subplots(figsize=(15, 15))
lnd.plotSites(fig, ax)
lnd.plotMigrationNetwork(fig, ax)
lnd.plotTraps(fig, ax)
lnd.plotTrapsNetwork(fig, ax)
ax.set_aspect('equal')