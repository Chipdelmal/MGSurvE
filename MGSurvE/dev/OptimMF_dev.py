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
import warnings
warnings.filterwarnings('ignore', 'The iteration is not making good progress')


(GENS, VERBOSE) = (1000, True)
if srv.isNotebook():
    (OUT_PTH, LND_TYPE, ID) = (
        '/home/chipdelmal/Documents/WorkSims/MGSurvE_Benchmarks/Sex/', 
        'UNIF', 'SX01'
)
else:
    (OUT_PTH, LND_TYPE, ID) = (argv[1], argv[2], argv[3].zfill(3))
###############################################################################
# Generating Pointsets
###############################################################################
if LND_TYPE == 'UNIF':
    ptsNum = 200
    bbox = ((-100, 100), (-80, 80))
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
###############################################################################
# Defining Traps
###############################################################################
traps = pd.DataFrame({
    'x': [0, 0],
    'y': [0, 0],
    't': [0, 0],
    'f': [0, 0]
})
tKernels = {
    'Male': {
        0: {'kernel': srv.exponentialDecay, 'params': {'A': .3, 'b': .1}}
    },
    'Female': {
        0: {'kernel': srv.exponentialDecay, 'params': {'A': .25, 'b': .075}}
    }
}
###############################################################################
# Defining Movement
###############################################################################
movementKernel = {
    'Male': {
        'kernelFunction': srv.zeroInflatedExponentialKernel,
        'kernelParams': {'params': [.075, 1.0e-10, math.inf], 'zeroInflation': .5}
    },
    'Female': {
        'kernelFunction': srv.zeroInflatedExponentialKernel,
        'kernelParams': {'params': [.075, 1.0e-10, math.inf], 'zeroInflation': .75}
    }
}
###############################################################################
# Setting Landscape Up
###############################################################################
lndM = srv.Landscape(
    points, traps=traps,
    kernelFunction=movementKernel['Male']['kernelFunction'],
    kernelParams=movementKernel['Male']['kernelParams'],
    trapsKernels=tKernels['Male']
)
lndF = srv.Landscape(
    points, traps=traps,
    kernelFunction=movementKernel['Female']['kernelFunction'],
    kernelParams=movementKernel['Female']['kernelParams'],
    trapsKernels=tKernels['Female']
)
srv.dumpLandscape(lndM, OUT_PTH, '{}_{}_M_CLN'.format(LND_TYPE, ID))
srv.dumpLandscape(lndF, OUT_PTH, '{}_{}_F_CLN'.format(LND_TYPE, ID))
# Needed auxiliary variables --------------------------------------------------
bbox = lndM.getBoundingBox()
trpMsk = srv.genFixedTrapsMask(lndM.trapsFixed)
###############################################################################
# Plot Landscape
###############################################################################
(fig, ax) = plt.subplots(1, 1, figsize=(15, 15), sharey=False)
lndM.plotSites(fig, ax, size=100)
lndM.plotMigrationNetwork(fig, ax, alphaMin=.3, lineWidth=50, lineColor='#03045e')
lndF.plotMigrationNetwork(fig, ax, alphaMin=.8, lineWidth=35, lineColor='#000000')
srv.plotClean(fig, ax, frame=False, bbox=bbox)
fig.savefig(
    path.join(OUT_PTH, '{}_{}_CLN.png'.format(LND_TYPE, ID)), 
    facecolor='w', bbox_inches='tight', pad_inches=0, dpi=300
)
###############################################################################
# GA Settings
############################################################################### 
POP_SIZE = int(10*(lndM.trapsNumber*1.25))
(MAT, MUT, SEL) = (
    {'mate': .3, 'cxpb': 0.5}, 
    {'mean': 0, 'sd': min([i[1]-i[0] for i in bbox])/5, 'mutpb': .4, 'ipb': .5},
    {'tSize': 3}
)
lndM_GA = deepcopy(lndM)
lndF_GA = deepcopy(lndF)
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
    trapsCoords=lndM_GA.trapsCoords, 
    fixedTrapsMask=trpMsk, coordsRange=bbox
)
toolbox.register("individualCreator", tools.initIterate, 
    creator.Individual, toolbox.initChromosome
)
toolbox.register("populationCreator", tools.initRepeat, 
    list, toolbox.individualCreator
)
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
    srv.calcSexFitness, 
    landscapeMale=lndM_GA,landscapeFemale=lndM_GA,
    maleWeight=.8, femaleWeight=1,
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
# Update with best results ----------------------------------------------------
minFits= logbook.select("min")
lndM.updateTrapsCoords(np.reshape(hof[0], (-1, 2)))
lndF.updateTrapsCoords(np.reshape(hof[0], (-1, 2)))
srv.dumpLandscape(lndM, OUT_PTH, '{}_{}_M_TRP'.format(LND_TYPE, ID))
srv.dumpLandscape(lndF, OUT_PTH, '{}_{}_F_TRP'.format(LND_TYPE, ID))
dta = pd.DataFrame(logbook)
srv.exportLog(logbook, OUT_PTH, '{}_{}_LOG'.format(LND_TYPE, ID))
###############################################################################
# Plot traps
############################################################################### 
(fig, ax) = plt.subplots(1, 1, figsize=(15, 15), sharey=False)
lndM.plotTraps(fig, ax, colors={0: '#f7258515'})
lndF.plotTraps(fig, ax, colors={0: '#32437615'})
srv.plotClean(fig, ax, frame=False, bbox=bbox)
srv.plotFitness(fig, ax, min(minFits))
fig.savefig(
    path.join(OUT_PTH, '{}_{}_TRP.png'.format(LND_TYPE, ID)), 
    facecolor='w', bbox_inches='tight', pad_inches=0, dpi=300
)
plt.close('all')
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
