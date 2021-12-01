#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
from copy import deepcopy
import matplotlib.pyplot as plt
from deap import base, creator, algorithms, tools
import MGSurvE as srv
from compress_pickle import dump, load


OUT_PTH = './'
###############################################################################
# Defining Landscape and Traps
###############################################################################
# pts = (
#     (0, 0, 0), 
#     (5, 10, 0), 
#     (10, 5, 0),
# )
# points = pd.DataFrame(pts, columns=('x', 'y', 't'))
ptsNum = 375
bbox = ((-100, 100), (-70, 70))
# xy = srv.ptsRegularGrid(ptsNum, bbox).T
xy = srv.ptsRandUniform(ptsNum, bbox).T
points = pd.DataFrame({'x': xy[0], 'y': xy[1], 't': [0]*xy.shape[1]})
# Traps info ------------------------------------------------------------------
traps = pd.DataFrame({
    'x': [0, 0, 0, 0, 0, 0], 
    'y': [0, 0, 0, 0, 0, 0], 
    't': [0, 0, 2, 1, 1, 1],
    'f': [0, 0, 0, 0, 0, 0]
})
tKernels = {
    0: {'kernel': srv.exponentialDecay, 'params': {'A': .75, 'b': .150}},
    1: {'kernel': srv.exponentialDecay, 'params': {'A': .50, 'b': .075}},
    2: {'kernel': srv.exponentialDecay, 'params': {'A': .25, 'b': .030}} 
}
###############################################################################
# Defining Landscape and Traps
###############################################################################
lnd = srv.Landscape(points, traps=traps, trapsKernels=tKernels)
srv.dumpLandscape(lnd, OUT_PTH, 'LND_ORG')
# lnd.calcFundamentalMatrix()
# lnd.getDaysTillTrapped()
bbox = lnd.getBoundingBox()
trpMsk = srv.genFixedTrapsMask(lnd.trapsFixed)
###############################################################################
# GA Settings
############################################################################### 
POP_SIZE = int(10*(lnd.trapsNumber*1.25))
(GENS, MAT, MUT, SEL) = (
    2000,
    {'mate': .5, 'cxpb': 0.5}, 
    {'mean': 0, 'sd': max([i[1]-i[0] for i in bbox])/4, 'mutpb': .2, 'ipb': .2},
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
srv.dumpLandscape(lnd, OUT_PTH, 'LND_OPT')
###############################################################################
# Plot
############################################################################### 
(fig, ax) = plt.subplots(1, 1, figsize=(15, 15), sharey=False)
lnd.plotSites(fig, ax, size=100)
lnd.plotMigrationNetwork(fig, ax)
lnd.plotTraps(fig, ax)
# lnd.plotTrapsNetwork(fig, ax)
srv.plotClean(fig, ax, frame=False)
ax.text(
    0.5, 0.5, '{:.3f}'.format(min(minFits)),
    horizontalalignment='center', verticalalignment='center',
    fontsize=100, color='#00000011',
    transform=ax.transAxes, zorder=5
)
fig.savefig(
    './GA2.png', facecolor='w',
    bbox_inches='tight', pad_inches=0, dpi=300
)
plt.close('all')