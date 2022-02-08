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
(ptsNum, trpsNum, bbox) = (200, 7, ((-100, 100), (-80, 80)))
gens = 250
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
            'params': [.075, 1.0e-10, math.inf], 'zeroInflation': .75
        }
    },
    'Female': {
        'kernelFunction': srv.zeroInflatedExponentialKernel,
        'kernelParams': {
            'params': [.075, 1.0e-10, math.inf], 'zeroInflation': .75
        }
    }
}
mMasks = {
    # Aquatic, Blood-Haunt, Sugar
    'Male': [
        [0.05, 0.35, 0.60],
        [0.45, 0.10, 0.45],
        [0.70, 0.10, 0.20],
    ],
    'Female': [
        [0.05, 0.70, 0.25],
        [0.30, 0.10, 0.60],
        [0.70, 0.10, 0.20],
    ]
}
###############################################################################
# Defining Traps
###############################################################################
nullTraps = [0]*trpsNum
traps = pd.DataFrame({
    'x': nullTraps, 'y': nullTraps, 
    'f': [True, False, False, False, False, False, False], 
    't': [0, 0, 0, 1, 1, 1, 1],
})
# 0: Ovitrap
# 1: CO2
tKernels = {
    'Male': {
        0: {'kernel': srv.exponentialDecay, 'params': {'A': 0.20, 'b': 0.250}},
        1: {'kernel': srv.exponentialDecay, 'params': {'A': 0.50, 'b': 0.075}}
    },
    'Female': {
        0: {'kernel': srv.exponentialDecay, 'params': {'A': 0.80, 'b': 0.115}},
        1: {'kernel': srv.exponentialDecay, 'params': {'A': 0.75, 'b': 0.075}}
    }
}
tMsk = {
    # Aquatic, Blood-Haunt, Sugar
    'Male': np.asarray([
        [1.0, 0.75, 0.50],
        [0.0, 1.0, 0.0]
    ]),
    'Female': np.asarray([
        [1.00, 0.50, 0.20],
        [0.50, 1.00, 0.20]
    ])
}
###############################################################################
# Setting Landscape Up
###############################################################################
lndM = srv.Landscape(
    points, traps=traps,
    maskingMatrix=mMasks['Male'],
    kernelFunction=movementKernel['Male']['kernelFunction'],
    kernelParams=movementKernel['Male']['kernelParams'],
    trapsKernels=tKernels['Male'], trapsMask=tMsk['Male'],
    trapsRadii=[.175, ]
)
lndF = srv.Landscape(
    points, traps=traps,
    maskingMatrix=mMasks['Female'],
    kernelFunction=movementKernel['Female']['kernelFunction'],
    kernelParams=movementKernel['Female']['kernelParams'],
    trapsKernels=tKernels['Female'], trapsMask=tMsk['Female'], 
    trapsRadii=[.175, ]
)
# Plotting male landscape .....................................................
(fig, ax) = plt.subplots(1, 1, figsize=(15, 15), sharey=False)
lndM.plotSites(fig, ax)
lndM.plotMigrationNetwork(fig, ax, alphaMin=.5, lineWidth=100)
srv.plotClean(fig, ax, frame=True)
fig.savefig(
    path.join(OUT_PTH, '{}_LND_M.png'.format(ID)), 
    facecolor='w', bbox_inches='tight', pad_inches=cst.pad, dpi=cst.dpi
)
plt.close('all')
# Plotting male landscape .....................................................
(fig, ax) = plt.subplots(1, 1, figsize=(15, 15), sharey=False)
lndF.plotSites(fig, ax)
lndF.plotMigrationNetwork(fig, ax, alphaMin=.5, lineWidth=100)
srv.plotClean(fig, ax, frame=True)
fig.savefig(
    path.join(OUT_PTH, '{}_LND_F.png'.format(ID)), 
    facecolor='w', bbox_inches='tight', pad_inches=cst.pad, dpi=cst.dpi
)
plt.close('all')
###############################################################################
# GA Settings
############################################################################### 
(weightMale, weightFemale) = (0.25, 0.75)
POP_SIZE = int(10*(lndM.trapsNumber*1.2))
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
# Plot GA
############################################################################### 
(fig, ax) = plt.subplots(figsize=(15, 15))
(fig, ax) = srv.plotGAEvolution(fig, ax, dta)
# srv.plotClean(fig, ax)
pthSave = path.join(OUT_PTH, '{}_GAP'.format(ID))
fig.savefig(
    pthSave,
    facecolor='w', bbox_inches='tight', pad_inches=.1, dpi=cst.dpi
)
###############################################################################
# Plot traps
############################################################################### 
(fig, ax) = plt.subplots(1, 1, figsize=(15, 15), sharey=False)
lndM.plotSites(fig, ax, size=150)
# Plot Networks ---------------------------------------------------------------
lndM.plotMigrationNetwork(fig, ax, alphaMin=.25, lineWidth=75, lineColor='#03045e')
lndF.plotMigrationNetwork(fig, ax, alphaMin=.25, lineWidth=75, lineColor='#03045e')
lndM.plotTraps(fig, ax, colors={0: '#3a0ca322', 1: '#a2d2ff33'}, lws=(2, 0), fill=True, ls=':', zorder=(10, 10))
lndF.plotTraps(fig, ax, colors={0: '#F966A922', 1: '#f1c0e833'}, lws=(2, 0), fill=True, ls='--', zorder=(25, 4))
# Other Stuff -----------------------------------------------------------------
srv.plotFitness(fig, ax, min(minFits), zorder=30)
srv.plotClean(fig, ax, frame=False, labels=False)
fig.savefig(
    path.join(OUT_PTH, '{}_TRP.png'.format(ID)), 
    facecolor='w', bbox_inches='tight', pad_inches=0.05, dpi=400
)
plt.close('all')
