#!/usr/bin/env python
# -*- coding: utf-8 -*-

import warnings
import numpy as np
import pandas as pd
from os import path
from sys import argv
from copy import deepcopy
import matplotlib.pyplot as plt
from compress_pickle import dump, load
from deap import base, creator, algorithms, tools
import Constants as cst
import MGSurvE as srv
warnings.filterwarnings('ignore', 'The iteration is not making good progress')


(GENS, VERBOSE, OUT_PTH) = (cst.gens, cst.verbose, cst.out_pth)
if srv.isNotebook():
    ID = 'Grid_LND_HOM'
else:
    ID = argv[1]
###############################################################################
# Load Landscape
###############################################################################
lnd = srv.loadLandscape(OUT_PTH, ID)
# Needed auxiliary variables --------------------------------------------------
(bbox, trpMsk) = (lnd.getBoundingBox(), srv.genFixedTrapsMask(lnd.trapsFixed))
###############################################################################
# GA Settings
############################################################################### 
TRPS_NUM = lnd.trapsCoords.shape[0]
POP_SIZE = int(10*(lnd.trapsNumber*1.25))
(MAT, MUT, SEL) = cst.gaParams
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
toolbox.register("initChromosome", srv.initDiscreteChromosome, 
    ptsIds=lndGA.pointID, 
    fixedTraps=lndGA.trapsFixed, 
    trapsSiteID=lndGA.trapsSiteID,
    banSites=lndGA.pointsTrapBanned
)
toolbox.register("individualCreator", tools.initIterate, 
    creator.Individual, toolbox.initChromosome
)
toolbox.register("populationCreator", tools.initRepeat, 
    list, toolbox.individualCreator
)
# Mate and mutate -------------------------------------------------------------
toolbox.register(
    "mate", srv.cxDiscreteUniform, 
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
    optimFunction=srv.getDaysTillTrappedPseudoInverse,
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
    stats=stats, halloffame=hof, verbose=True
)
# Update with best results ----------------------------------------------------
minFits= logbook.select("min")
lnd.updateTrapsCoords(np.reshape(hof[0], (-1, 2)))
srv.dumpLandscape(lnd, OUT_PTH, '{}_TRP'.format(ID))
dta = pd.DataFrame(logbook)
srv.exportLog(logbook, OUT_PTH, '{}_LOG'.format(ID))
###############################################################################
# Plot GA
############################################################################### 
(fig, ax) = plt.subplots(figsize=(15, 15))
(fig, ax) = srv.plotGAEvolution(fig, ax, dta)
# srv.plotClean(fig, ax)
pthSave = path.join(OUT_PTH, '{}_GAP-DOC'.format(ID))
fig.savefig(
    pthSave,
    facecolor='w', bbox_inches='tight', pad_inches=.1, dpi=cst.dpi
)
# Export plots ----------------------------------------------------------------
bbox = lnd.getBoundingBox()
trpMsk = srv.genFixedTrapsMask(lnd.trapsFixed)
(fig, ax) = plt.subplots(1, 1, figsize=(15, 15), sharey=False)
lnd.plotSites(fig, ax, size=200)
lnd.plotMaskedMigrationNetwork(fig, ax, alphaMin=.5, lineWidth=50)
lnd.plotTraps(fig, ax, size=200)
srv.plotClean(fig, ax, bbox=bbox, frame=False, pad=cst.pad_i)
srv.plotFitness(fig, ax, min(minFits), zorder=30)
fig.savefig(
    path.join(OUT_PTH, '{}_TRP-DOC.png'.format(ID)), 
    facecolor='w', bbox_inches='tight', pad_inches=cst.pad, dpi=cst.dpi
)
plt.close('all')