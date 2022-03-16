#!/usr/bin/env python
# -*- coding: utf-8 -*-

from sys import argv
import math
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
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import warnings
warnings.filterwarnings('ignore', 'The iteration is not making good progress')

(FXD_TRPS, TRPS_NUM, GENS) = (False, 10, 2000)
DIAG_VAL = 0.25
###############################################################################
# Debugging fixed traps at land masses
###############################################################################
OUT_PTH = './sims_out/'
if FXD_TRPS:
    ID = 'STP_FXD'
else:
    ID = 'STP_FXN'
###############################################################################
# Load Pointset
###############################################################################
SAO_TOME_LL = pd.read_csv(path.join('./GEO', 'STP_LatLon.csv'))
SAO_TOME_LL['t'] = [0]*SAO_TOME_LL.shape[0]
SAO_bbox = (
    (min(SAO_TOME_LL['x']), max(SAO_TOME_LL['x'])),
    (min(SAO_TOME_LL['y']), max(SAO_TOME_LL['y']))
)
SAO_cntr = [i[0]+(i[1]-i[0])/2 for i in SAO_bbox]
SAO_LIMITS = ((6.41, 6.79), (-0.0475, .45))
# Get location of minor land-masses -------------------------------------------
SAO_FIXED = [tuple(SAO_TOME_LL.loc[i][['x', 'y']]) for i in (24, 212)]
FXD_NUM = len(SAO_FIXED)
###############################################################################
# Load Migration Matrix
###############################################################################
migration = np.genfromtxt(
    path.join('./GEO', 'STP_Migration.csv'), delimiter=','
)
np.fill_diagonal(migration, DIAG_VAL)
SAO_TOME_MIG = normalize(migration, axis=1, norm='l1')
###############################################################################
# Defining Traps
#   North and South land masses (mini islands): 51, 239 (zero-indexed)
###############################################################################
(initTyp, initFxd) = ([0]*TRPS_NUM, [0]*TRPS_NUM)
(initLon, initLat) = ([SAO_cntr[0]]*TRPS_NUM, [SAO_cntr[1]]*TRPS_NUM)
if FXD_TRPS:
    for i in range(FXD_NUM):
        initFxd[TRPS_NUM-(i+1)] = 1
        initLon[TRPS_NUM-(i+1)] = SAO_FIXED[i][0]
        initLat[TRPS_NUM-(i+1)] = SAO_FIXED[i][1]
traps = pd.DataFrame({'x': initLon, 'y': initLat, 't': initTyp, 'f': initFxd})
tKer = {0: {'kernel': srv.exponentialDecay, 'params': {'A': .5, 'b': 100}}}
###############################################################################
# Setting Landscape Up
###############################################################################
lnd = srv.Landscape(
    SAO_TOME_LL, migrationMatrix=SAO_TOME_MIG,
    traps=traps, trapsKernels=tKer, landLimits=SAO_LIMITS
)
bbox = lnd.getBoundingBox()
trpMsk = srv.genFixedTrapsMask(lnd.trapsFixed)
###############################################################################
# Plot Results
###############################################################################
(fig, ax) = (
    plt.figure(figsize=(15, 15)),
    plt.axes(projection=ccrs.PlateCarree())
)
lnd.plotSites(fig, ax, size=250)
lnd.plotMigrationNetwork(
    fig, ax, lineWidth=10, alphaMin=.1, alphaAmplitude=2.5
)
lnd.plotLandBoundary(fig, ax)
srv.plotClean(fig, ax, bbox=lnd.landLimits)
fig.savefig(
    path.join(OUT_PTH, '{}_{:02d}_CLN.png'.format(ID, TRPS_NUM)), 
    facecolor='w', bbox_inches='tight', pad_inches=0.1, dpi=400
)
plt.close('all')
###############################################################################
# GA Settings
############################################################################### 
POP_SIZE = int(10*(lnd.trapsNumber*1.25))
(MAT, MUT, SEL) = (
    {'mate': .35, 'cxpb': 0.5}, 
    {
        'mean': 0, 
        'sd': max([abs(i[1]-i[0]) for i in bbox])/5, 
        'mutpb': .35, 'indpb': .5
    },
    {'tSize': 5}
)
VERBOSE = True
lndGA = deepcopy(lnd)
# Reducing the bbox for init sampling -----------------------------------------
redFract = 0
reduction = [(i[1]-i[0])/2*redFract for i in bbox]
bboxRed = [(i[0]+r, i[1]-r) for (i, r) in zip(bbox,reduction)]
###############################################################################
# Registering GA functions
############################################################################### 
toolbox = base.Toolbox()
creator.create(
    "FitnessMin", base.Fitness, 
    weights=(-1.0, )
)
creator.create(
    "Individual", list, 
    fitness=creator.FitnessMin
)
toolbox.register(
    "initChromosome", srv.initChromosome, 
    trapsCoords=lndGA.trapsCoords, 
    fixedTrapsMask=trpMsk, coordsRange=bboxRed
)
toolbox.register(
    "individualCreator", tools.initIterate, 
    creator.Individual, toolbox.initChromosome
)
toolbox.register(
    "populationCreator", tools.initRepeat, 
    list, toolbox.individualCreator
)
# Mate and mutate -------------------------------------------------------------
if FXD_TRPS:
    toolbox.register(
        "mutate", srv.mutateChromosome, 
        fixedTrapsMask=trpMsk, indpb=MUT['indpb'],
        randArgs={'loc': MUT['mean'], 'scale': MUT['sd']}
    )
    toolbox.register(
        "mate", srv.cxBlend, 
        fixedTrapsMask=trpMsk, alpha=MAT['mate']
    )
else:
    toolbox.register(
        "mutate", tools.mutGaussian, 
        mu=MUT['mean'], sigma=MUT['sd'], indpb=MUT['indpb']
    )
    toolbox.register(
        "mate", tools.cxBlend, 
        alpha=MAT['mate']
    )
# Select and evaluate ---------------------------------------------------------
toolbox.register(
    "select", tools.selTournament, 
    tournsize=SEL['tSize']
)
toolbox.register(
    "evaluate", srv.calcFitness, 
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
    pop, toolbox, 
    cxpb=MAT['cxpb'], mutpb=MUT['mutpb'], 
    ngen=GENS, stats=stats, halloffame=hof, verbose=VERBOSE
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
(fig, ax) = (
    plt.figure(figsize=(15, 15)),
    plt.axes(projection=ccrs.PlateCarree())
)
lnd.plotSites(fig, ax, size=250)
lnd.plotMigrationNetwork(
    fig, ax, lineWidth=10, alphaMin=.1, alphaAmplitude=2.5
)
lnd.plotTraps(fig, ax, zorders=(25, 20))
srv.plotFitness(fig, ax, min(dta['min']), fmt='{:.2f}')
lnd.plotLandBoundary(fig, ax)
srv.plotClean(fig, ax, bbox=lnd.landLimits)
fig.savefig(
    path.join(OUT_PTH, '{}_{:02d}_TRP.png'.format(ID, TRPS_NUM)), 
    facecolor='w', bbox_inches='tight', pad_inches=0.1, dpi=400
)
plt.close('all')


