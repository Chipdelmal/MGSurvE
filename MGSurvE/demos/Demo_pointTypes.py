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

(ID, OUT_PTH) = ('GA_DEMO_PT', './scratch/')
(PTS_NUM, PTS_TYPE, LND) = (150, 3, 'HOM')
###############################################################################
# Generate Sites
###############################################################################
bside = math.sqrt(PTS_NUM)*11.25
bbox = ((-bside, bside), (-bside, bside))
xy = srv.ptsRegularGrid(int(math.sqrt(PTS_NUM)), bbox).T
if LND == 'HOM':
    pType = [0]*xy.shape[1]
else:
    pType = np.random.choice(PTS_TYPE, xy.shape[1])
points = pd.DataFrame({'x': xy[0], 'y': xy[1], 't': pType})
###############################################################################
# Setup Points Movement
###############################################################################
msk = [
    [0.100, 0.700, 0.200],
    [0.100, 0.100, 0.800],
    [0.750, 0.125, 0.125],
]
movKer = {'params': [.075, 1.0e-10, math.inf], 'zeroInflation': .75}
###############################################################################
# Generate Traps
###############################################################################
trp = [
    [0, 0, 0, 0],
    [0, 0, 0, 0],
    [0, 0, 0, 0],
]
traps = pd.DataFrame(trp, columns=['x', 'y', 't', 'f'])
tker = {0: {'kernel': srv.exponentialDecay, 'params': {'A': 0.5, 'b': .05}}}
###############################################################################
# Setup Landscapes
###############################################################################
if LND == 'HOM':
    lnd = srv.Landscape(
        points, traps=traps, trapsKernels=tker
    )
else:
    lnd = srv.Landscape(
        points, maskingMatrix=msk, traps=traps, trapsKernels=tker
    )
bbox = lnd.getBoundingBox()
trpMsk = srv.genFixedTrapsMask(lnd.trapsFixed)
###############################################################################
# Plot Landscapes
###############################################################################
(fig, ax) = plt.subplots(1, 1, figsize=(10, 10), sharey=False)
lnd.plotSites(fig, ax, size=100)
lnd.plotMaskedMigrationNetwork(fig, ax, alphaMin=.6, lineWidth=25)
srv.plotClean(fig, ax, frame=False)
fig.savefig(
    path.join(OUT_PTH, '{}_{}.png'.format(ID, LND)), 
    facecolor='w', bbox_inches='tight', pad_inches=0.1, dpi=300
)
plt.close('all')
###############################################################################
# GA
###############################################################################
POP_SIZE = int(10*(lnd.trapsNumber*1.25))
(GENS, MAT, MUT, SEL) = (
    500,
    {'mate': .3, 'cxpb': 0.5}, 
    {'mean': 0, 'sd': min([i[1]-i[0] for i in bbox])/5, 'mutpb': .4, 'ipb': .5},
    {'tSize': 3}
)
VERBOSE = True
lndGA = deepcopy(lnd)
###############################################################################
# Registering GA functions
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
# Mate and mutate -------------------------------------------------------------
toolbox.register("mate", srv.cxBlend, 
    fixedTrapsMask=trpMsk, alpha=MAT['mate']
)
toolbox.register("mutate", srv.mutateChromosome, 
    fixedTrapsMask=trpMsk, 
    randArgs={'loc': MUT['mean'], 'scale': MUT['sd']}
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
bestTraps = np.reshape(hof[0], (-1, 2))
lnd.updateTrapsCoords(bestTraps)
srv.dumpLandscape(lnd, OUT_PTH, '{}_{}_TRP'.format(ID, LND))
dta = pd.DataFrame(logbook)
srv.exportLog(logbook, OUT_PTH, '{}_{}_LOG'.format(ID, LND))
###############################################################################
# Plot Landscape
############################################################################### 
(fig, ax) = plt.subplots(1, 1, figsize=(15, 15), sharey=False)
lnd.plotSites(fig, ax, size=100)
lnd.plotMaskedMigrationNetwork(fig, ax, alphaMin=.6, lineWidth=25)
lnd.plotTraps(fig, ax)
srv.plotClean(fig, ax, frame=False)
srv.plotFitness(fig, ax, min(dta['min']))
fig.savefig(
    path.join(OUT_PTH, '{}_{}_TRP'.format(ID, LND)), 
    facecolor='w', bbox_inches='tight', pad_inches=0.1, dpi=300
)
plt.close('all')