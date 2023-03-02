
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sympy.solvers import solve, nsolve
from sympy import Symbol
import MGSurvE as srv

LND = 3
ptsNum = 100
bbox = ((-100, 100), (-100, 100))
radii = (20, 100)
###############################################################################
# Select landscape type
###############################################################################
if LND == 0:
    xy = srv.ptsRegularGrid(ptsNum, bbox).T
elif LND == 1:
    xy = srv.ptsDonut(ptsNum, radii).T
elif LND == 2:
    xy = srv.ptsRandUniform(ptsNum, bbox).T
elif LND == 3:
    xy = srv.ptsRegularCircle(ptsNum, radii[-1]).T
###############################################################################
# Generate landscape object
###############################################################################
points = pd.DataFrame({'x': xy[0], 'y': xy[1], 't': [0]*xy.shape[1]})
lnd = srv.Landscape(points)
###############################################################################
# Plot
###############################################################################
(fig, ax) = plt.subplots(1, 1, figsize=(15, 15), sharey=False)
lnd.plotSites(fig, ax)
lnd.plotMigrationNetwork(fig, ax)
srv.plotClean(fig, ax, frame=False)