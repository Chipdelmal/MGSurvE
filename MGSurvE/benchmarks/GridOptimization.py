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
from timeit import default_timer as timer
import MGSurvE as srv
import warnings
warnings.filterwarnings('ignore', 'The iteration is not making good progress')


(TRP_NUM, PTS_NUM, REP) = (2, 100, 1)
PTH_O = '/home/chipdelmal/Documents/WorkSims/MGSurvE_Benchmarks'
GENS = 20
###############################################################################
# IDs and Timers
###############################################################################
expID = 'GB_{:02d}_{:05d}_{:02d}'.format(TRP_NUM, PTS_NUM, REP)
time = timer()
timers = {'start': time}
###############################################################################
# Generate Sites
###############################################################################
bside = math.sqrt(PTS_NUM)*11.25
bbox = ((-bside, bside), (-bside, bside))
xy = srv.ptsRegularGrid(int(math.sqrt(PTS_NUM)), bbox).T
points = pd.DataFrame({'x': xy[0], 'y': xy[1], 't': [0]*xy.shape[1]})
movKer = {'params': [.075, 1.0e-10, math.inf], 'zeroInflation': .75}
###############################################################################
# Generate Traps
###############################################################################
vLst = [0] * TRP_NUM
tKernels = {0: {'kernel': srv.exponentialDecay, 'params': {'A': .3, 'b': .05}}}
traps = pd.DataFrame({'x': vLst, 'y': vLst, 't': vLst, 'f': vLst})
###############################################################################
# Create Landscape
###############################################################################
lnd = srv.Landscape(
    points, kernelParams=movKer,
    traps=traps, trapsKernels=tKernels
)
srv.dumpLandscape(lnd, PTH_O, '{}_CLN'.format(expID))
lndGA = deepcopy(lnd)
###############################################################################
# GA Settings
############################################################################### 
POP_SIZE = int(10*(lnd.trapsNumber*1.25))
(MAT, MUT, SEL) = (
    {'mate': .3, 'cxpb': 0.5}, 
    {'mean': 0, 'sd': bside/5, 'mutpb': .5, 'ipb': .5},
    {'tSize': 3}
)
VERBOSE = False
# Timing ----------------------------------------------------------------------
timers['setup'] = timer()-timers['start']
###############################################################################
# Plot Landscape
###############################################################################
(fig, ax) = plt.subplots(1, 1, figsize=(15, 15), sharey=False)
lnd.plotSites(fig, ax, size=100)
lnd.plotMigrationNetwork(fig, ax, alphaMin=.6, lineWidth=50)
srv.plotClean(fig, ax, frame=False, bbox=bbox)
srv.saveFig(fig, ax, PTH_O, expID+'_CLN')
# Timing ----------------------------------------------------------------------
timers['plt_init'] = (timer()-timers['start'])+(timers['setup'])
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
    fixedTrapsMask=srv.genFixedTrapsMask(lnd.trapsFixed), coordsRange=bbox
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
minFits = logbook.select("min")
# Timing ----------------------------------------------------------------------
timers['ga'] = (timer()-timers['start'])+(timers['setup']+timers['plt_init'])
###############################################################################
# Plot final
############################################################################### 
lnd.plotTraps(fig, ax)
srv.plotClean(fig, ax, frame=False, bbox=bbox)
srv.plotFitness(fig, ax, min(minFits))
srv.saveFig(fig, ax, PTH_O, expID+'_TRP')
plt.close('all')
# Timing ----------------------------------------------------------------------
timers['plt_trp'] = (timer()-timers['start'])+(timers['setup']+timers['plt_init'])