#!/usr/bin/env python
# -*- coding: utf-8 -*-


from os import path
from glob import glob
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import MGSurvE as srv
plt.rc('font', size=15)
plt.rc('axes', titlesize=20)

(ID, PTH_O) = ('Grid', '/home/chipdelmal/Documents/WorkSims/MGSurvE_Benchmarks')
###############################################################################
# Load Experiments Data
###############################################################################
filePaths = glob(path.join(PTH_O, ID)+'/*TIM.csv')
fileNames = [path.basename(i) for i in filePaths]
fileExps = [tuple([int(j) for j in i.split('_')[1:3]]) for i in fileNames]
uniqueExps = list(set(fileExps))
###############################################################################
# Load Experiments Data
###############################################################################
exp = uniqueExps[0]
results = {}
for exp in uniqueExps:
    expPat = 'GB_{:02d}_{:05d}*TIM*'.format(*exp)
    expPaths = glob(path.join(path.join(PTH_O, ID), expPat))
    dta = [pd.read_csv(i).iloc[-1]['time']/60 for i in expPaths]
    summary = np.mean(dta)
    results[exp] = summary
###############################################################################
# Color palette
###############################################################################
cmap = srv.colorPaletteFromHexList([
    '#3a86ff', '#9fa0ff', '#9d4edd', '#000814'
])
yRange = (0, 60)
###############################################################################
# Scaling on Points
###############################################################################
(fig, ax) = plt.subplots(figsize=(10, 10))
points = sorted(set([i[0] for i in uniqueExps]))
traps = sorted(set([i[1] for i in uniqueExps]))
cols = [cmap(i) for i in np.arange(0, 1, 1/(len(traps)))]
for (ix, trap) in enumerate(traps):
    depen = [results[(i, trap)] for i in points]
    ax.plot(points, depen, label=trap, color=cols[ix])
    ax.set_ylim(*yRange)
    ax.set_xlim(min(points), max(points))
    plt.legend(title='Traps')
    plt.xlabel("Number of Sites")
    plt.ylabel("Time (minutes)")
    srv.saveFig(fig, ax, path.join(PTH_O, ID), 'PointsVTime')
###############################################################################
# Scaling on Traps
###############################################################################
(fig, ax) = plt.subplots(figsize=(10, 10))
cols = [cmap(i) for i in np.arange(0, 1, 1/(len(points)))]
for (ix, point) in enumerate(points):
    depen = [results[(point, i)] for i in traps]
    ax.plot(traps, depen, label=point, color=cols[ix])
    ax.set_ylim(*yRange)
    ax.set_xlim(min(traps), max(traps))
    plt.legend(title='Sites')
    plt.xlabel("Number of Traps")
    plt.ylabel("Time (minutes)")
    srv.saveFig(fig, ax, path.join(PTH_O, ID), 'TrapsVTime')