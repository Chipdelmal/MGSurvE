
import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from timeit import timeit
import MGSurvE as srv


(ID, TYPE, OUT_PTH) = ('PSO', 'Uniform', './demos_out/')
srv.makeFolder(OUT_PTH)

ptsNum = 25
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
mKer = {'params': srv.MEDIUM_MOV_EXP_PARAMS, 'zeroInflation': .25}
###############################################################################
# Traps
############################################################################### 
nullTraps = [0, 0, 0, 0]
traps = pd.DataFrame({
    'x': [0, 0, 0, 0], 
    'y': [0, 0, 0, 0], #[0, 0, 87.5, -87.5],
    't': [0, 1, 0, 1], 
    'f': [0, 0, 0, 0]
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
# Plot Kernel
###############################################################################
x = np.array(list(np.arange(1, 25, 1)))
(fig, ax) = plt.subplots(1, 1, figsize=(15, 5), sharey=False)
(fig, ax) = srv.plotMovementKernel(fig, ax, x, lnd)

