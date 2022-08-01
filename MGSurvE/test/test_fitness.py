
import math
import unittest
import numpy as np
import pandas as pd
import MGSurvE as srv

###############################################################################
# Setup Landscape
###############################################################################
ptsNum = 2
bbox = ((-20, 20), (-20, 20))
xy = srv.ptsRegularGrid(ptsNum, bbox).T
points = pd.DataFrame({
    'x': xy[0], 'y': xy[1], 
    't': [0]*xy.shape[1], 'id': range(0, xy.shape[1])
})
# Traps -----------------------------------------------------------------------
trapsNum = 2
nullTrap = [0]*trapsNum
traps = pd.DataFrame({
    'sid': nullTrap, 'x': [20, -20], 'y': [20, -20], 
    't': nullTrap, 'f': nullTrap
})
tKernels = {
    0: {'kernel': srv.exponentialDecay, 'params': {'A': .75, 'b': 0.1}}
}
# Instantiate Landscape -------------------------------------------------------
lnd = srv.Landscape(
    points, 
    kernelParams={'params': [.075, 1.0e-10, math.inf], 'zeroInflation': .75},
    traps=traps, trapsKernels=tKernels, pointsTrapBanned={5}, landLimits=bbox
)

# Traps at (0, 0)
lnd.updateTrapsCoords(np.asarray([[0, 0], [0, 0]]))
lnd.distanceMatrix
lnd.trapsDistances
lnd.migrationMatrix
lnd.maskedMigration
lnd.trapsMigration
# Traps at sites
# xy = srv.chromosomeIDtoXY((0, 3), lnd.pointID, lnd.pointCoords)
lnd.updateTrapsCoords(np.asarray([[-20, -20], [20, -20]]))
lnd.distanceMatrix
lnd.trapsDistances
lnd.migrationMatrix
lnd.trapsMigration