#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math
import pandas as pd
from sys import argv
import numpy as np
from os import path
from copy import deepcopy
import matplotlib.pyplot as plt
from deap import base, creator, algorithms, tools
import MGSurvE as srv
import warnings
warnings.filterwarnings('ignore', 'The iteration is not making good progress')

if srv.isNotebook():
    (OUT_PTH, LND_TYPE, ID) = ('./Lands', 'GRID', 'G03')
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
    ptsNum = 10
    bbox = ((-225, 225), (-225, 225))
    xy = srv.ptsRegularGrid(ptsNum, bbox).T
elif LND_TYPE == 'DNUT':
    ptsNum = 150
    radii = (100, 150)
    xy = srv.ptsDonut(ptsNum, radii).T
points = pd.DataFrame({
    'x': xy[0], 'y': xy[1], 
    't': [0]*xy.shape[1], 'id': range(0, xy.shape[1])
})
# Traps info ------------------------------------------------------------------
trapsNum = 8
nullTrap = [0]*trapsNum
tTypes = nullTrap[:]
tTypes[-1] = 1
traps = pd.DataFrame({
    'sid': nullTrap,
    'x': nullTrap, 'y': nullTrap,
    't': tTypes, 'f': nullTrap
})
tKernels = {
    0: {'kernel': srv.exponentialDecay, 'params': {'A': .75, 'b': 0.2}},
    1: {'kernel': srv.exponentialDecay, 'params': {'A': .75, 'b': 0.05}}
}
###############################################################################
# Defining Landscape and Traps
###############################################################################
lnd = srv.Landscape(
    points, kernelParams={'params': [.075, 1.0e-10, math.inf], 'zeroInflation': .75},
    traps=traps, trapsKernels=tKernels, pointsTrapBanned={5}
)
bbox = lnd.getBoundingBox()
trpMsk = srv.genFixedTrapsMask(lnd.trapsFixed)
###############################################################################
# GA Settings
############################################################################### 
POP_SIZE = int(10*(lnd.trapsNumber*1.25))
(GENS, MAT, MUT, SEL) = (
    500,
    {'cxpb':  0.3, 'indpb': 0.5}, 
    {'mutpb': 0.5, 'indpb': 0.35},
    {'tSize': 3}
)
VERBOSE = True
lndGA = deepcopy(lnd)
###############################################################################
# Implementing extension
###############################################################################
toolbox = base.Toolbox()
creator.create("FitnessMin", 
    base.Fitness, weights=(-1.0, )
)
creator.create("Individual", 
    list, fitness=creator.FitnessMin
)
# Init Population -------------------------------------------------------------
toolbox.register("initChromosome", srv.initDiscreteChromosome, 
    ptsIds=lndGA.pointID, 
    fixedTraps=lndGA.trapsFixed, 
    banSites=lndGA.pointsTrapBanned
)
toolbox.register("individualCreator", tools.initIterate, 
    creator.Individual, toolbox.initChromosome
)
toolbox.register("populationCreator", tools.initRepeat, 
    list, toolbox.individualCreator
)
# Mate and mutate -------------------------------------------------------------
toolbox.register("mate", srv.cxDiscreteUniform, 
    fixedTraps=lndGA.trapsFixed,
    indpb=MAT['indpb']
)
toolbox.register(
    "mutate", srv.mutateDiscreteChromosome,
    ptsIds=lndGA.pointID, 
    fixedTraps=lndGA.trapsFixed,
    indpb=MUT['indpb'],
    banSites=lndGA.pointsTrapBanned
)
# Select and evaluate ---------------------------------------------------------
toolbox.register("select", 
    tools.selTournament, tournsize=SEL['tSize']
)
toolbox.register("evaluate", 
    srv.calcDiscreteFitness, 
    landscape=lndGA,
    optimFunction=srv.getDaysTillTrapped,
    optimFunctionArgs={'outer': np.mean, 'inner': np.max}
)
###############################################################################
# Registering GA stats
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
bestChromosome = hof[0]
# Code this into a function ---------------------------------------------------
ptsIds = lndGA.pointID
siteIndex = [ptsIds.index(i) for i in bestChromosome]
trapXY = np.asarray([lndGA.pointCoords[i] for i in siteIndex])
# -----------------------------------------------------------------------------
lnd.updateTrapsCoords(trapXY)
dta = pd.DataFrame(logbook)
###############################################################################
# Plot Landscape
############################################################################### 
(fig, ax) = plt.subplots(1, 1, figsize=(15, 15), sharey=False)
lnd.plotSites(fig, ax, size=100)
lnd.plotMigrationNetwork(fig, ax, alphaMin=.6, lineWidth=25)
lnd.plotTraps(fig, ax)
srv.plotClean(fig, ax, frame=False)
srv.plotFitness(fig, ax, min(dta['min']))
fig.savefig(
    path.join(OUT_PTH, '{}_{}_TRP.png'.format(ID, LND_TYPE)), 
    facecolor='w', bbox_inches='tight', pad_inches=0.1, dpi=300
)
plt.close('all')


###############################################################################
# Drafts
###############################################################################
# trpsIDPos  = [0, 10, 55, 25]
# fixedTraps = [0, 0, 0, 1]
# trapsNum = lnd.trapsNumber
# ptsNum = lnd.pointNumber
# ptsIds = tuple((range(ptsNum)))

# chromB = srv.initDiscreteChromosome(lnd.pointID, lnd.trapsFixed, lnd.pointsTrapBanned)
# chromA = srv.mutateDiscreteChromosome(
#     chromB.copy(), lnd.pointID, lnd.trapsFixed, indpb=1
# )[0]
# print(chromA, chromB)
# print(srv.cxDiscreteUniform(chromA, chromB,  lnd.trapsFixed, indpb=.5))
# srv.calcDiscreteFitness(chromA, lnd)
# srv.calcDiscreteFitnessPseudoInverse(chromB, lnd)
# srv.calcDiscreteSexFitness(chromA, lnd, lnd)

# vct = [0]*100
# chrom = srv.initDiscreteChromosome(range(10), vct, {5, 6})
# set(chrom)

# (ub, chromSize) = (100, 100)
# chromA = srv.mutateDiscreteChromosome(
#     [0]*chromSize, range(1, ub), [0]*chromSize, indpb=1
# )[0]
# chromB = srv.mutateDiscreteChromosome(
#     [0]*chromSize, range(1, ub), [0]*chromSize, indpb=0
# )[0]
# noZero = (len([i for i in chromA if i==0]) == 0)
# allZero = (len([i for i in chromB if i==0]) == len(chromB))

# srv.cxDiscreteUniform(chromA, chromB, fixedTraps)