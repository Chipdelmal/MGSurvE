#!/usr/bin/env python
# -*- coding: utf-8 -*-

CORES = 16
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
import numpy as np
import pandas as pd
from os import path
from tqdm import tqdm
from sys import argv
from glob import glob
from random import shuffle
from random import randint
from time import perf_counter
from itertools import product
from SALib.sample import latin
import matplotlib.pyplot as plt
import matplotlib.ticker as tkr
from compress_pickle import dump, load
from termcolor import cprint
import auxiliary as aux
import MGSurvE as srv


if srv.isNotebook():
    DISCRETE = False
else:
    DISCRETE = int(argv[1])
# Constants for output --------------------------------------------------------
PTH_O = './sims_out/'
BBOX = ((-100, 100), (-100, 100))
srv.makeFolder(PTH_O)
# Experiment constants --------------------------------------------------------
(SEED, CORNERS) = (randint(0, 9999), True)
(GENS, REPS, LHS) = (500, 5, False)
(PTS_RAN, TRP_RAN, LAT_EXP) = ((50, 400, 50), (5, 30, 5), 100)
(BAN_PTS, BAN_TRP) = (None, None) # ((0, 200), (0, 20))
(SUM_STAT, INTERP) = (np.mean, 'linear')
###############################################################################
# Check for seed repetition
###############################################################################
app = ('DSC' if DISCRETE else 'CNT')
FILES = glob(path.join(PTH_O, f"timings_{app}*.bz"))
usedSeeds = set([int(i.split('-')[-1].split('.')[0]) for i in FILES])
while SEED in usedSeeds:
    SEED = randint(0, 9999)
###############################################################################
# Generate factorial tuples
###############################################################################
FACTORIAL = list(product(*[list(range(*PTS_RAN)), list(range(*TRP_RAN))]))
FACTORIAL = sorted(FACTORIAL, key=lambda x: x[1])
TIME = {}
###############################################################################
# Generate LHS tuples
###############################################################################
problem = {
    'num_vars': 2,
    'names': ['sites', 'traps'],
    'bounds': [[PTS_RAN[0], PTS_RAN[1]], [TRP_RAN[0], TRP_RAN[1]]]
}
LATIN = np.around(latin.sample(problem, LAT_EXP, seed=SEED)).astype(int)
if CORNERS:
    CORNERS = [
        [PTS_RAN[0], TRP_RAN[0]], [PTS_RAN[1], TRP_RAN[0]], 
        [PTS_RAN[0], TRP_RAN[1]], [PTS_RAN[1], TRP_RAN[1]], 
    ]
    LATIN = np.vstack((CORNERS, LATIN))
###############################################################################
# Iteration cycle
###############################################################################
# ix = -1
# (ptsNum, trpNum) = FACTORIAL[ix]
SAMPLE = (LATIN if LHS else FACTORIAL)
SAMPLE = FACTORIAL+[list(i) for i in LATIN]
# Filter ban range ------------------------------------------------------------
BAN_ENABLED = (BAN_PTS is not None) and (BAN_TRP is not None)
if BAN_ENABLED:
    BAN_RANGE = [
        (ptsNum in range(*BAN_PTS)) and (trpNum in range(*BAN_TRP)) 
        for (ptsNum, trpNum) in SAMPLE
    ]
    SAMPLE = [(s[0], s[1]) for (ix, s) in enumerate(SAMPLE) if not BAN_RANGE[ix]]
# Filter already run ----------------------------------------------------------
FILES = glob(path.join(PTH_O, f"timings_{app}*.bz"))
TIMES_LIST = [load(f) for f in FILES]
TIME = {k: v for d in TIMES_LIST for k, v in d.items()}
SAMPLE = [s for s in SAMPLE if not (tuple(s) in set(TIME))]
shuffle(SAMPLE)
# Iteration cycle -------------------------------------------------------------
cprint(
    f"* Running {len(SAMPLE)} experiments with {REPS} repetitions and {GENS} generations each on {app} setting!", 
    "red", "on_black"
)
for (ptsNum, trpNum) in tqdm(SAMPLE):
    # Setup sites -------------------------------------------------------------
    xy = srv.ptsRandUniform(ptsNum, BBOX).T
    pts = pd.DataFrame({'x': xy[0], 'y': xy[1], 't': [0]*xy.shape[1]})
    # Setup traps -------------------------------------------------------------
    nt= [0]*trpNum
    traps = pd.DataFrame({'sid': nt, 'x': nt, 'y': nt, 't': nt, 'f': nt})
    tKer = {0: {'kernel': srv.exponentialDecay, 'params': {'A': .75, 'b': .1}}}
    # Instantiate landscape ---------------------------------------------------
    lnd = srv.Landscape(pts, traps=traps, trapsKernels=tKer)
    # Time iterations ---------------------------------------------------------
    timings = []
    for rep in tqdm(range(REPS)):
        t0_rep = perf_counter()
        # Optimize ------------------------------------------------------------
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
    dump(TIME, path.join(PTH_O, f'timings_{app}-{SEED:04d}.bz2'))
###############################################################################
# Analyze resulting dictionary
###############################################################################
scale = 2
(XRAN, YRAN) = ((0, 30), (0, 400))
app = ('DSC' if DISCRETE else 'CNT')
title = ('discrete' if DISCRETE else 'continuous')
FILES = glob(path.join(PTH_O, f"timings_{app}*.bz"))
TIMES_LIST = [load(f) for f in FILES]
TIME = {k: v for d in TIMES_LIST for k, v in d.items()}
(XAXIS, YAXIS) = (
    {(i, 0): 0 for i in range(0, YRAN[1]+PTS_RAN[-1], int(PTS_RAN[-1]/1))},
    {(0, i): 0 for i in range(0, XRAN[1]+TRP_RAN[-1], int(TRP_RAN[-1]/1))}
)
TIME = TIME | XAXIS | YAXIS
(y, x) = np.array(list(TIME.keys())).T
z = np.array([scale*SUM_STAT(i)/60 for i in TIME.values()])
rs = aux.calcResponseSurface(x, y, z, mthd=INTERP)
(a, b) = ((min(x), max(x)), (min(y), max(y)))
(ran, rsG, rsS) = (rs['ranges'], rs['grid'], rs['surface'])
# Plot ------------------------------------------------------------------------
(cmin, cmax, cdelta) = (0, 90*scale, 10)
levels = np.arange(cmin, cmax+cdelta, cdelta)
if DISCRETE:
    cmap = srv.colorPaletteFromHexList(['#ffffff', '#8093f1', '#3a0ca3'])
else:
    cmap = srv.colorPaletteFromHexList(['#ffffff', '#8093f1', '#3a0ca3'])
    # ['#ffffff', '#ffafcc', '#ff006e'])
(lc, lw, ls) = ('#000000DD', 0.1, ":")
(fig, ax) = plt.subplots(figsize=(11, 10))
xy = ax.plot(rsG[0], rsG[1], 'k.', ms=1.25, alpha=.5, marker='o', zorder=10)
# cc = ax.contour(rsS[0], rsS[1], rsS[2], levels=levels, colors='#000000', linewidths=.5, alpha=1)
cs = ax.contourf(rsS[0], rsS[1], rsS[2], levels=levels, cmap=cmap, extend='max')
ax.set_xlabel("Number of Traps")
ax.set_ylabel("Number of Sites")
ax.set_title(f"Runtime over {scale*GENS} generations\n({title} optimization on {len(TIME)} samples)")
ax.vlines(list(set(x)), min(y), max(y), color=lc, lw=lw, ls=ls)
ax.hlines(list(set(y)), min(x), max(x), color=lc, lw=lw, ls=ls)
ax.vlines(range(5, max(x), 5), min(y), max(y), lw=lw*3, ls=ls)
ax.hlines(range(50, max(y), 50), min(x), max(x), lw=lw*3, ls=ls)
ax.set_xlim(*XRAN)
ax.set_ylim(*YRAN)
# ax.set_aspect('equal')
cbar = fig.colorbar(
    cs, ax=ax, ticks=np.linspace(0, max(z), 5), 
    format=tkr.FormatStrFormatter('%.f')
)
cbar.ax.set_ylabel('Time (minutes)')
cbar_ticks = np.linspace(cmin, cmax, num=len(levels), endpoint=True)
cbar.set_ticks(cbar_ticks)
fig.savefig(
    path.join(PTH_O, f'timings_{app}.png'), 
    facecolor='w', bbox_inches='tight', pad_inches=0.1, dpi=300
)