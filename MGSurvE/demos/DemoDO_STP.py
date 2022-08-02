#!/usr/bin/env python
# -*- coding: utf-8 -*-

from sys import argv
import numpy as np
import pandas as pd
from os import path
from sys import argv
from copy import deepcopy
import matplotlib.pyplot as plt
from numpy.random import uniform
from deap import base, creator, algorithms, tools
from sklearn.preprocessing import normalize
import MGSurvE as srv
import cartopy.crs as ccrs
import cartopy.feature as cfeature


if not srv.isNotebook():
    (FXD_TRPS, TRPS_NUM) = (int(argv[2]), int(argv[1]))
else:
    (FXD_TRPS, TRPS_NUM) = (False, 8)
###############################################################################
# Debugging fixed traps at land masses
###############################################################################
# OUT_PTH = '/Volumes/marshallShare/MGS_Benchmarks/STPVincenty/'
OUT_PTH = '/RAID5/marshallShare/MGS_Benchmarks/STPVincenty/'
# OUT_PTH = '/home/chipdelmal/Documents/WorkSims/MGSurvE_Benchmarks/STPVincenty'
if FXD_TRPS:
    ID = 'STP_DO_FXD'
else:
    ID = 'STP_DO_FXN'
GENS = 500
(IX_SPLIT, DIAG_VAL) = (27, 0.1)
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
SAO_cntr = [i[0]+(i[1]-i[0])/2 for i in SAO_bbox]
# SAO_TOME_LL = SAO_TOME_LL .rename(
#     columns={'lon': 'x', 'lat': 'y'}
# )
SAO_LIMITS = ((6.41, 6.79), (-0.0475, .45))
# Get location of minor land-masses -------------------------------------------
# SAO_FIXED = [tuple(SAO_TOME_LL.loc[i][['lon', 'lat']]) for i in (51, 239)]
SAO_FIXED = [51, 239]
FXD_NUM = len(SAO_FIXED)
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
#   North and South land masses (mini islands): 51, 239 (zero-indexed)
###############################################################################
(initTyp, initFxd) = ([0]*TRPS_NUM, [0]*TRPS_NUM)
(initLon, initLat) = ([
    uniform(*SAO_bbox[0], TRPS_NUM), uniform(*SAO_bbox[1], TRPS_NUM)
])
if FXD_TRPS:
    initFxd = ([0]*(TRPS_NUM-FXD_NUM) + [1]*FXD_NUM)
sid = [0]*(TRPS_NUM-FXD_NUM) + SAO_FIXED 
traps = pd.DataFrame({
    'sid': sid,
    'lon': initLon, 'lat': initLat, 
    't': initTyp, 'f': initFxd
})
tKer = {0: {'kernel': srv.exponentialDecay, 'params': {'A': 1, 'b': .0075}}}
###############################################################################
# Setting Landscape Up
###############################################################################
lnd = srv.Landscape(
    SAO_TOME_LL, migrationMatrix=SAO_TOME_MIG,
    traps=traps, trapsKernels=tKer, landLimits=SAO_LIMITS,
    trapsRadii=[.75, .5, .3],
)
bbox = lnd.getBoundingBox()
trpMsk = srv.genFixedTrapsMask(lnd.trapsFixed)
###############################################################################
# Plot Landscape
###############################################################################
# (fig, ax) = (
#     plt.figure(figsize=(15, 15)),
#     plt.axes(projection=ccrs.PlateCarree())
# )
# lnd.plotSites(fig, ax, size=250)
# # lnd.plotTraps(fig, ax)
# lnd.plotMigrationNetwork(
#     fig, ax, lineWidth=60, alphaMin=.1, alphaAmplitude=5,
# )
# lnd.plotLandBoundary(fig, ax)
# srv.plotClean(fig, ax, bbox=lnd.landLimits)
# fig.savefig(
#     path.join(OUT_PTH, '{}_{:02d}_CLN.png'.format(ID, TRPS_NUM)), 
#     facecolor='w', bbox_inches='tight', pad_inches=0.1, dpi=300
# )
# plt.close('all')
###############################################################################
# GA Settings
############################################################################### 
POP_SIZE = int(10*(lnd.trapsNumber*1.25))
(MAT, MUT, SEL) = (
    {'cxpb':  0.50, 'indpb': 0.35}, 
    {'mutpb': 0.45, 'indpb': 0.35},
    {'tSize': 3}
)
VERBOSE = True
lndGA = deepcopy(lnd)
# Reducing the bbox for init sampling -----------------------------------------
redFract = .25
reduction = [(i[1]-i[0])/2*redFract for i in bbox]
bboxRed = [(i[0]+r, i[1]-r) for (i, r) in zip(bbox,reduction)]
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
# Init Population -------------------------------------------------------------
toolbox.register("initChromosome", srv.initDiscreteChromosome, 
    ptsIds=lndGA.pointID, 
    fixedTraps=lndGA.trapsFixed, 
    # trapsSiteID=lndGA.trapsSiteID,
    banSites=lndGA.pointsTrapBanned
)
toolbox.register("individualCreator", tools.initIterate, 
    creator.Individual, toolbox.initChromosome
)
toolbox.register("populationCreator", tools.initRepeat, 
    list, toolbox.individualCreator
)
# Mate and mutate -------------------------------------------------------------
toolbox.register("mate", srv.cxDiscreteUniform, 
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
    optimFunctionArgs={'outer': np.mean, 'inner': np.mean}
)
###############################################################################
# Registering GA stats
############################################################################### 
pop = toolbox.populationCreator(n=POP_SIZE)
hof = tools.HallOfFame(1)
stats = tools.Statistics(lambda ind: ind.fitness.values)   
stats.register("min", np.min)
stats.register("median", np.median)
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
dta = pd.DataFrame(logbook)
srv.dumpLandscape(lnd, OUT_PTH, '{}_{:02d}_TRP'.format(ID, TRPS_NUM), fExt='pkl')
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
    fig, ax, lineWidth=60, alphaMin=.1, alphaAmplitude=5
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
