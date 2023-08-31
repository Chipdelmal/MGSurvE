
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
step = ('zeroInflation' in lnd.kernelParams.keys())
x = np.arange(0, 20, 1)
y = [0]*x.shape[0]
xy = np.array([x, y]).T
distMat = srv.calcDistanceMatrix(xy, distFun=math.dist)
kernMat = lnd.kernelFunction(distMat, **lnd.kernelParams)
(fig, ax) = plt.subplots(1, 1, figsize=(15, 5), sharey=False)
if not step:
    ax.plot(x, kernMat, lw=4)
    ax.scatter(x, kernMat[0], zorder=5, marker='x', s=50)
else:
    ax.plot(x[1:], kernMat[0][1:], lw=4, color='#8093f1')
    ax.plot([0, 0], [0, kernMat[0][0]], lw=8, color='#ec0868')
    ax.vlines(x[1:], ymin=0, ymax=1, zorder=5, lw=0.5, alpha=0.25, color='#7371fc')
ax.set_xlim(x[0], x[-1])
ax.set_ylim(0, 0.5) # kernMat[0, 0])