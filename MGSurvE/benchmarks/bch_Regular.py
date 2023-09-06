#!/usr/bin/env python
# -*- coding: utf-8 -*-

CORES = 8
###############################################################################
# Load libraries and limit cores
###############################################################################
import os
os.environ["OMP_NUM_THREADS"] = str(CORES)
os.environ["OPENBLAS_NUM_THREADS"] = str(CORES)
os.environ["MKL_NUM_THREADS"] = str(CORES)
os.environ["VECLIB_MAXIMUM_THREADS"] = str(CORES)
os.environ["NUMEXPR_NUM_THREADS"] = str(CORES)
# Load libraries --------------------------------------------------------------
import warnings
warnings.filterwarnings("ignore")
import warnings
import numpy as np
import pandas as pd
from os import path
from tqdm import tqdm
from time import perf_counter
from itertools import product
import matplotlib.pyplot as plt
from compress_pickle import dump, load
from termcolor import colored, cprint
import auxiliary as aux
import MGSurvE as srv

# Constants for output --------------------------------------------------------
PTH_O = './sims_out/'
BBOX = ((-100, 100), (-100, 100))
srv.makeFolder(PTH_O)
# Experiment constants --------------------------------------------------------
(GENS, REPS, DISCRETE) = (1000, 10, True)
(PTS_RAN, TRP_RAN) = ((1, 10, 2), (1, 10, 2))
SUM_STAT = np.median
###############################################################################
# Generate factorial tuples
###############################################################################
FACTORIAL = list(product(*[
    list(range(*PTS_RAN)), 
    list(range(*TRP_RAN))
]))
TIME = {}
###############################################################################
# Iteration cycle
###############################################################################
ix = 0
(ptsNum, trpNum) = FACTORIAL[ix]
cprint(f"* Running {len(FACTORIAL)} experiments. Please wait!", "red", "on_black")
for (ptsNum, trpNum) in tqdm(FACTORIAL):
    # Setup sites -------------------------------------------------------------
    xy = srv.ptsRandUniform(ptsNum, BBOX).T
    pts = pd.DataFrame({'x': xy[0], 'y': xy[1], 't': [0]*xy.shape[1]})
    # Setup traps -------------------------------------------------------------
    nt= [0]*trpNum
    traps = pd.DataFrame({'sid': nt, 'x': nt, 'y': nt, 't': nt, 'f': nt})
    tKer = {0: {'kernel': srv.exponentialDecay, 'params': {'A': .75, 'b': .1}}}
    # Instantiate landscape ---------------------------------------------------
    lnd = srv.Landscape(pts, traps=traps, trapsKernels=tKer)
    # Optimize ----------------------------------------------------------------
    timings = []
    for rep in range(REPS):
        t0_rep = perf_counter()
        if not DISCRETE:
            (lnd, logbook) = srv.optimizeTrapsGA(
                lnd, generations=GENS, pop_size='auto', verbose=False,
                mating_params='auto', mutation_params='auto', 
                selection_params='auto', fitFuns={'inner': np.sum, 'outer': np.max}
            )
        else:
            (lnd, logbook) = srv.optimizeDiscreteTrapsGA(
                lnd, generations=GENS, pop_size='auto', verbose=False,
                mating_params='auto', mutation_params='auto', 
                selection_params='auto', fitFuns={'inner': np.sum, 'outer': np.max}
            )
        tf_rep = perf_counter()
        timings.append(tf_rep-t0_rep)
    TIME[(lnd.pointNumber, lnd.trapsNumber)] = timings
    # Export results ----------------------------------------------------------
    app = ('DSC' if DISCRETE else 'CNT')
    dump(TIME, path.join(PTH_O, f'timings_{app}.bz2'))
###############################################################################
# Analyze resulting dictionary
###############################################################################
app = ('DSC' if DISCRETE else 'CNT')
TIME = load(path.join(PTH_O, f'timings_{app}.bz2'))
cmap = srv.colorPaletteFromHexList(['#ffffff', '#8093f1'])
# Plot ------------------------------------------------------------------------
(x, y) = np.array(list(TIME.keys())).T
z = np.array([SUM_STAT(i) for i in TIME.values()])
rs = aux.calcResponseSurface(x, y, z)
(a, b) = ((min(x), max(x)), (min(y), max(y)))
(ran, rsG, rsS) = (rs['ranges'], rs['grid'], rs['surface'])
(fig, ax) = plt.subplots(figsize=(10, 10))
xy = ax.plot(rsG[0], rsG[1], 'k.', ms=5, alpha=.75, marker='x')
cc = ax.contour(rsS[0], rsS[1], rsS[2], colors='#000000', linewidths=.5, alpha=1)
cs = ax.contourf(rsS[0], rsS[1], rsS[2], cmap=cmap, extend='max')
ax.set_xlabel("Number of Sites")
ax.set_ylabel("Number of Traps")
ax.set_title(f"Runtime over {GENS} generations ({app})")
ax.set_aspect('equal')
cbar = fig.colorbar(cs, ax=ax, ticks=np.linspace(0, 1, 5))
cbar.ax.set_ylabel('Time (minutes)')
fig.savefig(
    path.join(PTH_O, f'timings_{app}.png'), 
    facecolor='w', bbox_inches='tight', pad_inches=0, dpi=300
)