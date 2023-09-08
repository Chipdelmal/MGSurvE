
import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from timeit import timeit
import MGSurvE as srv


(ID, TYPE) = ('PSO', 'Uniform')
###############################################################################
# Setting
############################################################################### 
ptsNum = 200
radii = (75, 100)
pTypesProb =[0.05, 0.70, 0.25]
bbox = ((-200, 200), (-150, 150))
###############################################################################
# Pointset
############################################################################### 
if TYPE == 'Grid':
    (ptsNum, ptsTypes) = (int(math.sqrt(ptsNum)), len(pTypesProb))
    xy = srv.ptsRegularGrid(ptsNum, (bbox[0], bbox[0])).T
elif TYPE == 'Uniform':
    (ptsNum, ptsTypes) = (ptsNum, len(pTypesProb))
    xy = srv.ptsRandUniform(ptsNum, bbox).T
elif TYPE == 'Ring':
    (ptsNum, radii, ptsTypes) = (ptsNum, radii, len(pTypesProb))
    xy = srv.ptsDonut(ptsNum, radii).T
points = pd.DataFrame({'x': xy[0], 'y': xy[1], 't': [0]*xy.shape[1]})
mKer = {'params': [.075, 1.0e-10, math.inf], 'zeroInflation': .75}
###############################################################################
# Traps
############################################################################### 
nullTraps = [0, 0, 0, 0]
traps = pd.DataFrame({
    'x': [0, 0, 0, 0], 'y': [0, 0, 0, 0],
    't': [0, 1, 0, 1], 'f': [0, 0, 0, 0]
})
tKer = {
    0: {'kernel': srv.exponentialDecay, 'params': {'A': .75, 'b': .100}},
    1: {'kernel': srv.exponentialDecay, 'params': {'A': .50, 'b': .050}}
}
###############################################################################
# Landscape
############################################################################### 
lnd = srv.Landscape(points, kernelParams=mKer, traps=traps, trapsKernels=tKer)
###############################################################################
# Test and Compare
############################################################################### 
(tau, sitesN, trapsN, iters) = (
    lnd.trapsMigration, ptsNum, traps.shape[0], 50000
)
# Equivalency -----------------------------------------------------------------
mA = np.sum(srv.getFundamentalMatrix(tau, sitesN, trapsN), axis=1)
mB = np.sum(srv.getFundamentalMatrixPseudoInverse(tau, sitesN, trapsN), axis=1)
mC = srv.getFundamentalVector(tau, sitesN)
equivalency = all([
    all(np.isclose(mA, mB)), all(np.isclose(mA, mC)), all(np.isclose(mB, mC))
])
print("Methods are equivalent?: {}".format(equivalency))
# Timing ----------------------------------------------------------------------
tA = timeit(lambda: np.sum(srv.getFundamentalMatrix(tau, sitesN, trapsN), axis=1), number=iters)
tB = timeit(lambda: np.sum(srv.getFundamentalMatrixPseudoInverse(tau, sitesN, trapsN), axis=1), number=iters)
tC = timeit(lambda: srv.getFundamentalVector(tau, sitesN), number=iters)
print("* Time inverse: {}\n* Time pseudo: {}\n* Time solve: {}\nOver {} iterations!".format(tA, tB, tC, iters))