#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
from copy import deepcopy
import matplotlib.pyplot as plt
from deap import base, creator, algorithms, tools
import MGSurvE as srv


###############################################################################
# Defining Landscape and Traps
###############################################################################
pts = (
    (0, 0, 0), 
    (5, 10, 0), 
    (10, 5, 0),
)
points = pd.DataFrame(pts, columns=('x', 'y', 't'))
# Traps info ------------------------------------------------------------------
traps = pd.DataFrame({
    'x': [0.5, 2.0], 
    'y': [0.0, 2.0], 
    't': [0, 1],
    'f': [0, 0]
})
tKernels = {
    0: {'kernel': srv.exponentialDecay, 'params': {'A': .30, 'b': 2}},
    1: {'kernel': srv.exponentialDecay, 'params': {'A': .50, 'b': 1}} 
}
###############################################################################
# Defining Landscape and Traps
###############################################################################
lnd = srv.Landscape(points, traps=traps, trapsKernels=tKernels)
lnd.calcFundamentalMatrix()
lnd.getDaysTillTrapped()
bbox = lnd.getBoundingBox()
trpMsk = srv.genFixedTrapsMask(lnd.trapsFixed)
###############################################################################
# GA Settings
############################################################################### 
POP_SIZE = int(10*(lnd.trapsNumber*1.25))
(GENS, MAT, MUT, SEL) = (
    1000,
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
toolbox.register("initChromosome", 
    srv.initChromosome, trapsNum=lndGA.trapsNumber, coordsRange=bbox
)
toolbox.register("individualCreator", 
    tools.initIterate, creator.Individual, toolbox.initChromosome
)
toolbox.register("populationCreator", 
    tools.initRepeat, list, toolbox.individualCreator
)
# Custom ---------------------------------------------------------------------
toolbox.register("mate", 
    srv.cxBlend, fixedTrapsMask=trpMsk, 
    alpha=MAT['mate']
)
toolbox.register("mutate", 
    srv.mutateChromosome, fixedTrapsMask=trpMsk, 
    randArgs={'loc': MUT['mean'], 'scale': MUT['sd']}
)
# Original --------------------------------------------------------------------
# toolbox.register(
#     "mate", tools.cxBlend, 
#     alpha=MAT['mate']
# )
# toolbox.register(
#     "mutate", tools.mutGaussian, 
#     mu=MUT['mean'], sigma=MUT['sd'], indpb=MUT['ipb']
# )
# Select and evaluate ---------------------------------------------------------
toolbox.register("select", 
    tools.selTournament, tournsize=SEL['tSize']
)
toolbox.register("evaluate", 
    srv.calcFitness, 
    landscape=lndGA, trapsKernels=tKernels,
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

###############################################################################
# Dev
############################################################################### 
# srv.calcFitness(np.asarray(traps[['x', 'y']]), landscape=lnd, trapsKernels=tKernels)
# srv.calcFitness(np.asarray(traps[['x', 'y']]), landscape=lndGA, trapsKernels=tKernels)
# lndGA = deepcopy(lnd)
# traps = pd.DataFrame({
#     'x': [0.5, 3.0, 2.0], 
#     'y': [0.0, 0.0, 2.0], 
#     't': [0, 1, 0],
#     'f': [1, 1, 0]
# })
# srv.calcFitness(np.asarray(traps[['x', 'y']]), landscape=lndGA, trapsKernels=tKernels)
# srv.calcFitness(np.asarray(traps[['x', 'y']]), landscape=lnd, trapsKernels=tKernels)
# lndGA