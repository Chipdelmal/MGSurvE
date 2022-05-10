import operator
import random

from deap import base
from deap import creator
from deap import tools
import unittest

import math
import numpy as np
from os import path
import matplotlib.pyplot as plt
import pandas as pd
from copy import deepcopy
from deap import base, creator, algorithms, tools
import MGSurvE as srv

(ID, OUT_PTH) = ('PSO_testing', './test/')

def test_movable_pso():
     # defining landscape ----------------------------------------------------
    ptsNum = 100
    radii = (75, 100)
    xy = srv.ptsDonut(ptsNum, radii).T
    points = pd.DataFrame({'x': xy[0], 'y': xy[1], 't': [0]*xy.shape[1]})
    mKer = {'params': [.075, 1.0e-10, math.inf], 'zeroInflation': .75}
    # defining traps --------------------------------------------------------
    nullTraps = [0, 0, 0, 0]
    traps = pd.DataFrame({
        'x': nullTraps, 'y': nullTraps,
        't': nullTraps, 'f': nullTraps
    })
    tKer = {0: {'kernel': srv.exponentialDecay, 'params': {'A': .5, 'b': .1}}}
    # setting landscape up --------------------------------------------------------
    lnd = srv.Landscape(
        points, kernelParams=mKer,
        traps=traps, trapsKernels=tKer
    )
    bbox = lnd.getBoundingBox()
    trpMsk = srv.genFixedTrapsMask(lnd.trapsFixed)
    # (fig, ax) = plt.subplots(1, 1, figsize=(15, 15), sharey=False)
    # lnd.plotSites(fig, ax, size=100)
    # lnd.plotMigrationNetwork(fig, ax, alphaMin=.6, lineWidth=25)
    # lnd.plotTraps(fig, ax)
    # srv.plotClean(fig, ax, frame=False)
    # fig.savefig(
    #     path.join(OUT_PTH, '{}_TRP1.png'.format(ID)), 
    #     facecolor='w', bbox_inches='tight', pad_inches=0.1, dpi=300
    # )
    # GA settings -------------------------------------------------------------------
    POP_SIZE = int(10*(lnd.trapsNumber*1.25))
    (GENS, MAT, MUT, SEL) = (
        200,
        {'mate': .3, 'cxpb': 0.5}, 
        {'mean': 0, 'sd': min([i[1]-i[0] for i in bbox])/5, 'mutpb': .5, 'ipb': .5},
        {'tSize': 3}
    )
    VERBOSE = True
    lndGA = deepcopy(lnd)
    # PSO -------------------------------------------------------------------
    pso = srv.Particle_Swarm(traps,-150, 150, lnd, 5, 100)
    # self, traps, p_min, p_max, lnd, num_particles=50, num_gens=500, s_min=-3, s_max=3

    # (pop, logbook, best) = pso.evaluate()
    # Reshaping the chromosome and plotting ---------------------------------------------
    # bestTraps = np.reshape(best, (-1, 2))
    # lnd.updateTrapsCoords(bestTraps)
    # Generate the plot -----------------------------------------------------------
    # (fig, ax) = plt.subplots(1, 1, figsize=(15, 15), sharey=False)
    # lnd.plotSites(fig, ax, size=100)
    # lnd.plotMigrationNetwork(fig, ax, alphaMin=.6, lineWidth=25)
    # lnd.plotTraps(fig, ax)
    # srv.plotClean(fig, ax, frame=False, bbox=bbox)
    # fig.savefig(
    #     './images/pso.png',
    #     facecolor='w', bbox_inches='tight', 
    #     pad_inches=.1, dpi=300
    # )
   
test_movable_pso()

def test_immovable_pso():
    # defining landscape ----------------------------------------------------
    ptsNum = 100
    radii = (75, 100)
    xy = srv.ptsDonut(ptsNum, radii).T
    points = pd.DataFrame({'x': xy[0], 'y': xy[1], 't': [0]*xy.shape[1]})
    mKer = {'params': [.075, 1.0e-10, math.inf], 'zeroInflation': .75}
    
    # defining traps --------------------------------------------------------
    traps = pd.DataFrame({
    'x': [0, 0, 0, 0], 
    'y': [0, 0, 87.5, -87.5],
    't': [0, 1, 0, 1], 
    'f': [0, 0, 1, 1]})

    num_traps = traps.shape[0] 

    tKer = {
    0: {'kernel': srv.exponentialDecay, 'params': {'A': .5, 'b': .1}},
    1: {'kernel': srv.exponentialDecay, 'params': {'A': .25, 'b': .05}}
    }

    # setting landscape up --------------------------------------------------------
    lnd = srv.Landscape(
        points, kernelParams=mKer,
        traps=traps, trapsKernels=tKer
    )
    bbox = lnd.getBoundingBox()
    trpMsk = srv.genFixedTrapsMask(lnd.trapsFixed)

    # Initialize an instance of PSO
    opt = srv.Particle_Swarm(traps,-150, 150, lnd, 5, 50)

    # (fig, ax) = plt.subplots(1, 1, figsize=(15, 15), sharey=False)
    # lnd.plotSites(fig, ax, size=100)
    # lnd.plotMigrationNetwork(fig, ax, alphaMin=.6, lineWidth=25)
    # lnd.plotTraps(fig, ax)
    # srv.plotClean(fig, ax, frame=False)
    # fig.savefig(
    #     path.join(OUT_PTH, '{}_TRP1.png'.format(ID)), 
    #     facecolor='w', bbox_inches='tight', pad_inches=0.1, dpi=300
    # )

    # Initialize an instance of PSO
    opt = srv.Particle_Swarm(traps,-150, 150, lnd, 5, 50)
    # self, traps, p_min, p_max, lnd, num_particles=50, num_gens=500, s_min=-3, s_max=3

    # Run PSO 
    # (pop, logbook, best) = opt.evaluate() 

    ###############################################################################
    # Reshaping the chromosome and plotting 
    ###############################################################################
    # bestTraps = np.reshape(best, (-1, 2))
    # lnd.updateTrapsCoords(bestTraps)

    # Generate the plot -----------------------------------------------------------
    # (fig, ax) = plt.subplots(1, 1, figsize=(15, 15), sharey=False)
    # lnd.plotSites(fig, ax, size=100)
    # lnd.plotMigrationNetwork(fig, ax, alphaMin=.6, lineWidth=25)
    # lnd.plotTraps(fig, ax)
    # srv.plotClean(fig, ax, frame=False, bbox=bbox)
    # fig.savefig(
    #     './images/pso.png',
    #     facecolor='w', bbox_inches='tight', 
    #     pad_inches=.1, dpi=300
    # )

test_immovable_pso()

###############################################################################
# Main
###############################################################################
if __name__ == '__main__':
    unittest.main()

