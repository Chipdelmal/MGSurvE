#!/usr/bin/env python
# -*- coding: utf-8 -*-

from sys import argv
import math
from vincenty import vincenty
import numpy as np
import pandas as pd
from os import path
from sys import argv
from copy import deepcopy
import matplotlib.pyplot as plt
from deap import base, creator, algorithms, tools
from compress_pickle import dump, load
from sklearn.preprocessing import normalize
import MGSurvE as srv
import warnings
warnings.filterwarnings('ignore', 'The iteration is not making good progress')


(ID, OUT_PTH) = (
    'STP', 
    # '/RAID5/marshallShare/MGS_Benchmarks/STP/'
    '/home/chipdelmal/Documents/WorkSims/MGSurvE_Benchmarks/STPVincenty/'
)
TRPS_NUM = 3 # int(argv[1]) # 3
IX_SPLIT = 27
DIAG_VAL = 0
GENS = 200
###############################################################################
# Load Pointset
###############################################################################
sites = pd.read_csv(path.join(OUT_PTH, 'GEO', 'stp_cluster_sites_pop_v5_fixed.csv'))
sites['t'] = [0]*sites.shape[0]
SAO_TOME_LL = sites.iloc[IX_SPLIT:]
SAO_bbox = (
    (min(SAO_TOME_LL['lon']), max(SAO_TOME_LL['lon'])),
    (min(SAO_TOME_LL['lat']), max(SAO_TOME_LL['lat']))
)
SAO_TOME_LL = SAO_TOME_LL .rename(
    columns={'lon': 'x', 'lat': 'y'}
)
###############################################################################
# Load Migration Matrix
###############################################################################
migration = np.genfromtxt(
    path.join(OUT_PTH, 'GEO', 'kernel_cluster_v6a.csv'), delimiter=','
)
msplit = migration[IX_SPLIT:,IX_SPLIT:]
np.fill_diagonal(msplit, DIAG_VAL)
SAO_TOME_MIG = normalize(msplit, axis=1, norm='l1')
###############################################################################
# Defining Traps
###############################################################################
nullTraps = [0]*TRPS_NUM
(lonTrap, latTrap) = (
    np.random.uniform(SAO_bbox[0][0], SAO_bbox[0][1], TRPS_NUM),
    np.random.uniform(SAO_bbox[1][0], SAO_bbox[1][1], TRPS_NUM)
)
traps = pd.DataFrame({
    'x': lonTrap, 'y': latTrap,
    't': nullTraps, 'f': nullTraps
})
tKer = {0: {'kernel': srv.exponentialDecay, 'params': {'A': .5, 'b': 100}}}
###############################################################################
# Setting Landscape Up
###############################################################################
lnd = srv.Landscape(
    SAO_TOME_LL, migrationMatrix=SAO_TOME_MIG,
    traps=traps, trapsKernels=tKer,
    distanceFunction=math.dist
)
bbox = lnd.getBoundingBox()
trpMsk = srv.genFixedTrapsMask(lnd.trapsFixed)
###############################################################################
# Plot Landscape
###############################################################################
(fig, ax) = plt.subplots(1, 1, figsize=(15, 15), sharey=False)
lnd.plotSites(fig, ax)
# lnd.plotTraps(fig, ax)
lnd.plotMigrationNetwork(
    fig, ax, 
    lineWidth=5, alphaMin=.5, alphaAmplitude=2.5,
)
srv.plotClean(fig, ax, frame=False, labels=False)
fig.savefig(
    path.join(OUT_PTH, '{}_{:02d}_CLN.png'.format(ID, TRPS_NUM)), 
    facecolor='w', bbox_inches='tight', pad_inches=0.1, dpi=300
)
###############################################################################
# GA Settings
############################################################################### 
POP_SIZE = int(10*(lnd.trapsNumber*1.25))
(MAT, MUT, SEL) = (
    {'mate': .35, 'cxpb': 0.5}, 
    {
        'mean': 0, 
        'sd': min([abs(i[1]-i[0]) for i in bbox])/5, 
        'mutpb': .35, 'ipb': .5
    },
    {'tSize': 5}
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
bestTraps = np.reshape(bestChromosome, (-1, 2))
lnd.updateTrapsCoords(bestTraps)
srv.dumpLandscape(lnd, OUT_PTH, '{}_{:02d}_TRP'.format(ID, TRPS_NUM))
dta = pd.DataFrame(logbook)
srv.exportLog(logbook, OUT_PTH, '{}_{:02d}_LOG'.format(ID, TRPS_NUM))
###############################################################################
# Plot Results
###############################################################################
(fig, ax) = plt.subplots(1, 1, figsize=(15, 15), sharey=False)
lnd.plotSites(fig, ax)
lnd.plotMigrationNetwork(
    fig, ax, 
    lineWidth=5, alphaMin=.5, alphaAmplitude=5,
)
lnd.plotTraps(fig, ax, zorders=(25, 20))
srv.plotFitness(fig, ax, min(dta['min']), fmt='{:.5f}')
srv.plotClean(fig, ax, frame=False, labels=False)
fig.savefig(
    path.join(OUT_PTH, '{}_{:02d}_TRP.png'.format(ID, TRPS_NUM)), 
    facecolor='w', bbox_inches='tight', pad_inches=0.1, dpi=200
)
plt.close('all')
###############################################################################
# Debug
############################################################################### 
# bestChromosome = [6.65, 0]
# bestChromosome = srv.initChromosome(
#     lndGA.trapsCoords, fixedTrapsMask=trpMsk, coordsRange=bbox
# )
# bestTraps = np.reshape(bestChromosome, (-1, 2))
# fit = srv.calcFitness(
#     bestChromosome,
#     landscape=lnd,
#     optimFunction=srv.getDaysTillTrapped,
#     optimFunctionArgs={'outer': np.mean, 'inner': np.mean}
# )
# lnd.updateTrapsCoords(bestTraps)
# (fig, ax) = plt.subplots(1, 1, figsize=(15, 15), sharey=False)
# lnd.plotSites(fig, ax)
# lnd.plotTraps(fig, ax)
# srv.plotFitness(fig, ax, fit[0])
# srv.plotClean(fig, ax, frame=True, labels=True)
# fig.savefig(
#     path.join(OUT_PTH, '{}_DBG.png'.format(ID)), 
#     facecolor='w', bbox_inches='tight', pad_inches=0.1, dpi=300
# )
# srv.initChromosome(
#     lndGA.trapsCoords, fixedTrapsMask=trpMsk, coordsRange=bbox
# )
# bbox
