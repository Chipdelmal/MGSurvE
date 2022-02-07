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
import Constants as cst
import MGSurvE as srv
import warnings
warnings.filterwarnings('ignore', 'The iteration is not making good progress')


(OUT_PTH, ID) = (cst.out_pth, 'SX')
(ptsNum, trpsNum, bbox) = (200, 6, ((-100, 100), (-80, 80)))
gens = 1000
###############################################################################
# Generating Pointsets
###############################################################################
xy = srv.ptsRandUniform(ptsNum, bbox).T
points = pd.DataFrame({'x': xy[0], 'y': xy[1], 't': [0]*xy.shape[1]})
# Sites' Resources: {0: Aquatic, 1: Blood haunt, 2: Sugar source}
resourceProbs = [0.2, 0.5, 0.3]
points['t'] = np.random.choice(len(resourceProbs), ptsNum, p=resourceProbs)
###############################################################################
# Defining Movement
###############################################################################
movementKernel = {
    'Male': {
        'kernelFunction': srv.zeroInflatedExponentialKernel,
        'kernelParams': {
            'params': [.050, 1.0e-10, math.inf], 'zeroInflation': .5
        }
    },
    'Female': {
        'kernelFunction': srv.zeroInflatedExponentialKernel,
        'kernelParams': {
            'params': [.025, 1.0e-10, math.inf], 'zeroInflation': .7
        }
    }
}
###############################################################################
# Defining Traps
###############################################################################
nullTraps = [0]*trpsNum
traps = pd.DataFrame({
    'x': nullTraps, 'y': nullTraps, 
    'f': [False, False, False, True, False, False], 
    't': [0, 1, 0, 1, 0, 1],
})
tKernels = {
    'Male': {
        0: {'kernel': srv.exponentialDecay, 'params': {'A': .25, 'b': .125}},
        1: {'kernel': srv.exponentialDecay, 'params': {'A': .25, 'b': .125}}
    },
    'Female': {
        0: {'kernel': srv.exponentialDecay, 'params': {'A': .5, 'b': .100}},
        1: {'kernel': srv.exponentialDecay, 'params': {'A': .75, 'b': .075}}
    }
}
tMasks = {
    'Male': [
        [0.20, 0.10, 0.70],
        [0.40, 0.20, 0.40],
        [0.80, 0.05, 0.15]
    ],
    'Female': [
        [0.10, 0.70, 0.20],
        [0.30, 0.10, 0.60],
        [0.70, 0.10, 0.20]
    ]
}
###############################################################################
# Setting Landscape Up
###############################################################################
lndM = srv.Landscape(
    points, traps=traps,
    maskingMatrix=tMasks['Male'],
    kernelFunction=movementKernel['Male']['kernelFunction'],
    kernelParams=movementKernel['Male']['kernelParams'],
    trapsKernels=tKernels['Male'], trapsRadii=[.1, ]
)
lndF = srv.Landscape(
    points, traps=traps,
    maskingMatrix=tMasks['Female'],
    kernelFunction=movementKernel['Female']['kernelFunction'],
    kernelParams=movementKernel['Female']['kernelParams'],
    trapsKernels=tKernels['Female'], trapsRadii=[.1, ]
)
# Plotting male landscape .....................................................
(fig, ax) = plt.subplots(1, 1, figsize=(15, 15), sharey=False)
lndM.plotSites(fig, ax)
lndM.plotMigrationNetwork(fig, ax, alphaMin=.75, lineWidth=30)
srv.plotClean(fig, ax, frame=True)
fig.savefig(
    path.join(OUT_PTH, '{}_LND_M.png'.format(ID)), 
    facecolor='w', bbox_inches='tight', pad_inches=cst.pad, dpi=cst.dpi
)
plt.close('all')
# Plotting male landscape .....................................................
(fig, ax) = plt.subplots(1, 1, figsize=(15, 15), sharey=False)
lndF.plotSites(fig, ax)
lndF.plotMigrationNetwork(fig, ax, alphaMin=.75, lineWidth=90)
srv.plotClean(fig, ax, frame=True)
fig.savefig(
    path.join(OUT_PTH, '{}_LND_F.png'.format(ID)), 
    facecolor='w', bbox_inches='tight', pad_inches=cst.pad, dpi=cst.dpi
)
plt.close('all')
###############################################################################
# GA Settings
############################################################################### 
(weightMale, weightFemale) = (.5, 1)
POP_SIZE = int(10*(lndM.trapsNumber*1.25))
(GENS, MAT, MUT, SEL, VERBOSE) = (
    gens,
    {'mate': .3, 'cxpb': 0.5}, 
    {'mean': 0, 'sd': min([i[1]-i[0] for i in bbox])/5, 'mutpb': .4, 'ipb': .5},
    {'tSize': 3},
    True
)
lndM_GA = deepcopy(lndM)
lndF_GA = deepcopy(lndF)
# Needed auxiliary variables --------------------------------------------------
bbox = lndM.getBoundingBox()
trpMsk = srv.genFixedTrapsMask(lndM.trapsFixed)
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
    "mate", srv.cxBlend,
    fixedTrapsMask=trpMsk,
    alpha=MAT['mate']
)
toolbox.register(
    "mutate", srv.mutateChromosome, 
    fixedTrapsMask=trpMsk, 
    randArgs={'loc': MUT['mean'], 'scale': MUT['sd']}
)
# Select and evaluate ---------------------------------------------------------
toolbox.register("select", 
    tools.selTournament, tournsize=SEL['tSize']
)
toolbox.register("evaluate", 
    srv.calcSexFitness, 
    landscapeMale=lndM_GA,landscapeFemale=lndF_GA,
    weightMale=weightMale, weightFemale=weightFemale,
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
srv.exportLandscape(lndM, OUT_PTH, ID+'M')
srv.exportLandscape(lndF, OUT_PTH, ID+'F')
dta = pd.DataFrame(logbook)
###############################################################################
# Plot traps
############################################################################### 
(fig, ax) = plt.subplots(1, 1, figsize=(15, 15), sharey=False)
lndM.plotSites(fig, ax, size=100)
# Plot Networks ---------------------------------------------------------------
lndM.plotMigrationNetwork(fig, ax, alphaMin=.3, lineWidth=50, lineColor='#03045e')
lndF.plotMigrationNetwork(fig, ax, alphaMin=.3, lineWidth=35, lineColor='#03045e')
lndF.plotTraps(fig, ax, colors={0: '#f7258522', 1: '#f7258522'}, lws=(2, 0), fill=True, ls='--', zorder=(25, 4))
lndM.plotTraps(fig, ax, colors={0: '#a06cd522', 1: '#a06cd522'}, lws=(2, 0), fill=True, ls=':', zorder=(30, 5))
# Other Stuff -----------------------------------------------------------------
srv.plotFitness(fig, ax, min(minFits), zorder=30)
srv.plotClean(fig, ax, frame=False, bbox=bbox, labels=False)
fig.savefig(
    path.join(OUT_PTH, '{}_TRP.png'.format(ID)), 
    facecolor='w', bbox_inches='tight', pad_inches=0.05, dpi=300
)
plt.close('all')
