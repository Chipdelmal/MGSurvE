
import math
import numpy as np
from os import path
import matplotlib.pyplot as plt
import pandas as pd
import MGSurvE as srv


(ID, TYPE, OUT_PTH) = ('PSO', 'Ring', './demos_out/')
srv.makeFolder(OUT_PTH)


gens = 2000
ptsNum = 500
radii = (425, 500)
pTypesProb =[0.25, 0.75]
bbox = ((-500, 500), (-350, 350))
pad = 10
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
pType = np.random.choice(ptsTypes, xy.shape[1])
points = pd.DataFrame({'x': xy[0], 'y': xy[1], 't': pType}) # [0]*xy.shape[1]})
mKer = {'params': [.075, 1.0e-10, math.inf], 'zeroInflation': .75}
###############################################################################
# Traps
############################################################################### 
traps = pd.DataFrame({
    'x': [0, 0, 0, 0, 0, 0], 
    'y': [0, 0, 0, 0, 0, 0],
    't': [0, 1, 0, 1, 0, 1], 
    'f': [0, 0, 0, 0, 0, 0]
})
tKer = {
    0: {'kernel': srv.exponentialDecay, 'params': {'A': .75, 'b': .050}},
    1: {'kernel': srv.exponentialDecay, 'params': {'A': .50, 'b': .030}}
}
###############################################################################
# Landscape
############################################################################### 
lnd = srv.Landscape(
    points, kernelParams=mKer,
    traps=traps, trapsKernels=tKer
)
bbox = lnd.getBoundingBox()
trpMsk = srv.genFixedTrapsMask(lnd.trapsFixed)
###############################################################################
# Plot Landscape
############################################################################### 
(fig, ax) = plt.subplots(1, 1, figsize=(15, 15), sharey=False)
lnd.plotSites(fig, ax, size=100)
lnd.plotMigrationNetwork(fig, ax, alphaMin=.6, lineWidth=25)
srv.plotClean(fig, ax, frame=False, bbox=bbox, pad=(5, 5))
fig.savefig(
    path.join(OUT_PTH, '{}_{}-CLN.png'.format(ID, TYPE)),
    facecolor='w', bbox_inches='tight', 
    pad_inches=1, dpi=300
)
plt.close('all')
###############################################################################
# PSO
############################################################################### 
(GENS, PARTS, SPD, PHI) = (
    gens,
    traps.shape[0]*20,
    (-max(max(bbox))/30, max(max(bbox))/30),
    (max(max(bbox))/15, max(max(bbox))/15)
)
pso = srv.Particle_Swarm(
    lnd=lnd,
    traps=traps,
    num_particles=PARTS, num_gens=GENS, 
    p_min=min(bbox[0][0], bbox[1][0]), p_max=max(bbox[1][0], bbox[1][1]),  
    s_min=SPD[0], s_max=SPD[1],
    phi1=PHI[0], phi2=PHI[1],
    optimFunctionArgs={'outer': np.max, 'inner': np.sum}
)
(pop, logbook, _) = pso.evaluate()
best = list(logbook[logbook['min']==min(logbook['min'])]['traps'])[0]
bestTraps = np.reshape(best, (-1, 2))
lnd.updateTrapsCoords(bestTraps)
###############################################################################
# Export Results
############################################################################### 
dta = pd.DataFrame(logbook)
srv.dumpLandscape(lnd, OUT_PTH, '{}_{}-TRP'.format(ID, TYPE), fExt='pkl')
srv.exportLog(logbook, OUT_PTH, '{}_{}-LOG'.format(ID, TYPE))
###############################################################################
# Plot Results
############################################################################### 
(fig, ax) = plt.subplots(1, 1, figsize=(15, 15), sharey=False)
lnd.plotSites(fig, ax, size=100)
lnd.plotMigrationNetwork(fig, ax, alphaMin=.6, lineWidth=25)
lnd.plotTraps(fig, ax)
srv.plotFitness(fig, ax, min(logbook['min']), zorder=30)
srv.plotClean(fig, ax, frame=False, bbox=bbox, pad=(5, 5))
fig.savefig(
    path.join(OUT_PTH, '{}_{}-TRP.png'.format(ID, TYPE)),
    facecolor='w', bbox_inches='tight', 
    pad_inches=1, dpi=300
)
plt.close('all')
