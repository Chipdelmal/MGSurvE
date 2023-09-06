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
import math
import warnings
import numpy as np
import pandas as pd
from os import path
from sys import argv
from copy import deepcopy
from time import perf_counter
from itertools import product
import matplotlib.pyplot as plt
from compress_pickle import dump, load
import MGSurvE as srv

# Constants for output --------------------------------------------------------
PTH_O = './sims_out/'
BBOX = ((-100, 100), (-100, 100))
srv.makeFolder(PTH_O)
# Experiment constants --------------------------------------------------------
(GENS, REPS) = (10, 5)
(PTS_RAN, TRP_RAN) = ((3, 10, 2), (1, 5, 1))
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
for (ptsNum, trpNum) in FACTORIAL:
    # Setup sites -------------------------------------------------------------
    xy = srv.ptsRegularGrid(ptsNum, (BBOX[0], BBOX[0])).T
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
        (lnd, logbook) = srv.optimizeTrapsGA(
            lnd, generations=GENS, pop_size='auto', verbose=False,
            mating_params='auto', mutation_params='auto', 
            selection_params='auto', fitFuns={'inner': np.sum, 'outer': np.max}
        )
        tf_rep = perf_counter()
        timings.append(tf_rep-t0_rep)
    TIME[(ptsNum, trpNum)] = timings
    dump(TIME, path.join(PTH_O, 'timings.bz2'))
###############################################################################
# Analyze resulting dictionary
###############################################################################
load(path.join(PTH_O, 'timings.bz2'))