#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math
import numpy as np
import pandas as pd
from os import path
from sys import argv
from copy import deepcopy
import matplotlib.pyplot as plt
from deap import base, creator, algorithms, tools
from compress_pickle import dump, load
import MGSurvE as srv
import networkx as nx
from sklearn.preprocessing import normalize
import numpy as np
import warnings
warnings.filterwarnings('ignore', 'The iteration is not making good progress')
# scp -r lab:/RAID5/marshallShare/MGS_Demos/* '/home/chipdelmal/Documents/GitHub/MGSurvE/MGSurvE/dev/Lands'

if srv.isNotebook():
    (OUT_PTH, LND_TYPE, ID) = ('./Lands', 'DNUT', 'D01')
else:
    (OUT_PTH, LND_TYPE, ID) = (argv[1], argv[2], argv[3].zfill(3))
###############################################################################
# Defining Landscape and Traps
###############################################################################
if LND_TYPE == 'UNIF':
    ptsNum = 400
    bbox = ((-225, 225), (-175, 175))
    xy = srv.ptsRandUniform(ptsNum, bbox).T
elif LND_TYPE == 'GRID':
    ptsNum = 20
    bbox = ((-225, 225), (-225, 225))
    xy = srv.ptsRegularGrid(ptsNum, bbox).T
elif LND_TYPE == 'DNUT':
    ptsNum = 300
    radii = (100, 150)
    xy = srv.ptsDonut(ptsNum, radii).T
points = pd.DataFrame({'x': xy[0], 'y': xy[1], 't': [0]*xy.shape[1]})
# Traps info ------------------------------------------------------------------
traps = pd.DataFrame({
    'x': [0, 0, 0, 0, 0],
    'y': [0, 0, 0, 0, 0],
    't': [0, 0, 0, 0, 0],
    'f': [0, 0, 0, 0, 0]
})
tKernels = {
    0: {'kernel': srv.exponentialDecay, 'params': {'A': .3, 'b': .05}},
    1: {'kernel': srv.exponentialDecay, 'params': {'A': .35, 'b': .04}},
    2: {'kernel': srv.exponentialDecay, 'params': {'A': .25,  'b': .025}} ,
    3: {'kernel': srv.sigmoidDecay,     'params': {'A': .2, 'rate': .5, 'x0': 1}}
}
###############################################################################
# Defining Landscape and Traps
###############################################################################
lnd = srv.Landscape(
    points, kernelParams={'params': [.075, 1.0e-10, math.inf], 'zeroInflation': .75},
    traps=traps, trapsKernels=tKernels
)
srv.dumpLandscape(lnd, OUT_PTH, '{}_{}_CLN'.format(LND_TYPE, ID))
# lnd.calcFundamentalMatrix()
# lnd.getDaysTillTrapped()
bbox = lnd.getBoundingBox()
trpMsk = srv.genFixedTrapsMask(lnd.trapsFixed)
###############################################################################
# Plot Landscape
###############################################################################
(fig, ax) = plt.subplots(1, 1, figsize=(15, 15), sharey=False)
lnd.plotSites(fig, ax, size=100)
lnd.plotMigrationNetwork(fig, ax, alphaMin=.6, lineWidth=50)
srv.plotClean(fig, ax, frame=False, bbox=bbox)
fig.savefig(
    path.join(OUT_PTH, '{}_{}_CLN.png'.format(LND_TYPE, ID)), 
    facecolor='w', bbox_inches='tight', pad_inches=0, dpi=300
)
plt.close('all')
###############################################################################
# GA Settings
############################################################################### 
POP_SIZE = int(10*(lnd.trapsNumber*1.25))
(GENS, MAT, MUT, SEL) = (
    25,
    {'mate': .3, 'cxpb': 0.5}, 
    {'mean': 0, 'sd': min([i[1]-i[0] for i in bbox])/5, 'mutpb': .5, 'ipb': .5},
    {'tSize': 3}
)
VERBOSE = True
lndGA = deepcopy(lnd)
###############################################################################
# Registering Functions for GA
############################################################################### 
toolbox = base.Toolbox()
creator.create("FitnessMin", 
    base.Fitness, weights=(-1.0, )
)
creator.create("Individual", 
    list, fitness=creator.FitnessMin
)
toolbox.register("initChromosome", srv.initChromosome, 
    trapsCoords=lndGA.trapsCoords, 
    fixedTrapsMask=trpMsk, coordsRange=bbox
)
toolbox.register("individualCreator", tools.initIterate, 
    creator.Individual, toolbox.initChromosome
)
toolbox.register("populationCreator", tools.initRepeat, 
    list, toolbox.individualCreator
)
# Custom ---------------------------------------------------------------------
toolbox.register(
    "mate", tools.cxBlend, 
    alpha=MAT['mate']
)
toolbox.register(
    "mutate", tools.mutGaussian, 
    mu=MUT['mean'], sigma=MUT['sd'], indpb=MUT['ipb']
)
# Select and evaluate ---------------------------------------------------------
toolbox.register("select", 
    tools.selTournament, tournsize=SEL['tSize']
)
toolbox.register("evaluate", 
    srv.calcFitness, 
    landscape=lndGA,
    optimFunction=srv.getDaysTillTrapped,
    optimFunctionArgs={'outer': np.mean, 'inner': np.max}
)
###############################################################################
# Registering functions for GA stats
############################################################################### 
pop = toolbox.populationCreator(n=POP_SIZE)
hof = tools.HallOfFame(1)
stats = tools.Statistics(lambda ind: ind.fitness.values)   
stats.register("min", np.min)
stats.register("avg", np.mean)
stats.register("max", np.max)
stats.register("best", lambda fitnessValues: fitnessValues.index(min(fitnessValues)))
stats.register("traps", lambda fitnessValues: pop[fitnessValues.index(min(fitnessValues))])
###############################################################################
# Optimization Cycle
############################################################################### 
(pop, logbook) = algorithms.eaSimple(
    pop, toolbox, cxpb=MAT['cxpb'], mutpb=MUT['mutpb'], ngen=GENS, 
    stats=stats, halloffame=hof, verbose=VERBOSE
)
###############################################################################
# Get and Export Results
############################################################################### 
(maxFits, meanFits, bestIndx, minFits, traps) = logbook.select(
    "max", "avg", "best", "min", "traps"
)
lnd.updateTrapsCoords(np.reshape(hof[0], (-1, 2)))
srv.dumpLandscape(lnd, OUT_PTH, '{}_{}_TRP'.format(LND_TYPE, ID))
dta = pd.DataFrame(logbook)
srv.exportLog(logbook, OUT_PTH, '{}_{}_LOG'.format(LND_TYPE, ID))
###############################################################################
# Plot Landscape
############################################################################### 
(fig, ax) = plt.subplots(1, 1, figsize=(15, 15), sharey=False)
lnd.plotSites(fig, ax, size=100)
lnd.plotMigrationNetwork(fig, ax, alphaMin=.6, lineWidth=25)
lnd.plotTraps(fig, ax)
srv.plotClean(fig, ax, frame=False, bbox=bbox)
srv.plotFitness(fig, ax, min(minFits))
fig.savefig(
    path.join(OUT_PTH, '{}_{}_TRP.png'.format(LND_TYPE, ID)), 
    facecolor='w', bbox_inches='tight', pad_inches=0, dpi=300
)
plt.close('all')
print("Done!")
###############################################################################
# Plot GA
############################################################################### 
(fig, ax) = plt.subplots(figsize=(15, 15))
(fig, ax) = srv.plotGAEvolution(fig, ax, dta)
# srv.plotClean(fig, ax)
pthSave = path.join(
    OUT_PTH, '{}_{}_GAP'.format(LND_TYPE, ID)
)
fig.savefig(
    pthSave,
    facecolor='w', bbox_inches='tight', pad_inches=.1, dpi=300
)


# (trpsNum, dims) = (12, 2)
# trpsFxd = [0]*trpsNum
# initMsk = [True] * trpsNum * dims
# chrom = np.random.uniform(-10, 10, trpsNum*2)
# trpMsk = srv.genFixedTrapsMask(trpsFxd)

# org = np.copy(chrom)

# mod = srv.mutateChromosomeAsymmetric(
#     np.copy(chrom), trpMsk, 
#     randArgs={
#         'x': {'loc': 10, 'scale': 100}, 
#         'y': {'loc': 0, 'scale': 0}
#     }
# )
# totalX = np.sum([org[a]!=mod[0][a] for a in range(len(org))])

# mod = srv.mutateChromosomeAsymmetric(
#     np.copy(chrom), trpMsk, 
#     randArgs={
#         'x': {'loc': 0, 'scale': 0}, 
#         'y': {'loc': 10, 'scale': 100}
#     }
# )
# totalY = np.sum([org[a]!=mod[0][a] for a in range(len(org))])

# mod = srv.mutateChromosomeAsymmetric(
#     np.copy(chrom), trpMsk, 
#     randArgs={
#         'x': {'loc': 10, 'scale': 0}, 
#         'y': {'loc': 10, 'scale': 1}
#     }
# )
# totalXY = np.sum([org[a]!=mod[0][a] for a in range(len(org))])
