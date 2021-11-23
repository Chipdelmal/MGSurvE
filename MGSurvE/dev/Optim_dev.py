#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import MGSurvE as srv
import matplotlib.pyplot as plt

###############################################################################
# Defining Landscape and Traps
###############################################################################
pts = (
    (0.0, 0.0, 0), 
    (2.0, 0.5, 0), 
    (2.5, 1.5, 0),
)
points = pd.DataFrame(pts, columns=('x', 'y', 't'))
# Traps info ------------------------------------------------------------------
traps = pd.DataFrame({
    'x': [0.5, 3.0, 2.0], 
    'y': [0.0, 0.0, 2.0], 
    't': [0, 1, 0],
    'f': [1, 1, 0]
})
tKernels = {
    0: {'kernel': srv.exponentialDecay, 'params': {'A': .30, 'b': 2}},
    1: {'kernel': srv.exponentialDecay, 'params': {'A': .50, 'b': 1}} 
}
# Land creation ---------------------------------------------------------------
lnd = srv.Landscape(points, traps=traps, trapsKernels=tKernels)
lnd.calcFundamentalMatrix()
lnd.getDaysTillTrapped()
###############################################################################
# Defining Landscape and Traps
###############################################################################
chrom = srv.initChromosome(lnd.trapsNumber, coordsRange=((0, 5), (0, 5)))
fxdTrpsMsk = srv.genFixedTrapsMask(lnd.trapsFixed)
mut = srv.mutateChromosome(chrom, fxdTrpsMsk)

chrom, mut

(pre1, pre2) = (chrom.copy(), mut.copy())
(ind1, ind2) = srv.cxBlend(chrom, mut, fxdTrpsMsk)

[np.isclose(a, b) for (a, b) in zip(pre1, ind1)]
[np.isclose(a, b) for (a, b) in zip(pre2, ind2)]



###############################################################################
# Mutation function test
###############################################################################
trpsNum = 10
trpsFxd = [0]*trpsNum
dims = 2

initChrom = srv.initChromosome(trpsNum, coordsRange=((0, 5), (0, 5)))

results = []
for i in range(trpsNum):
    trpsFxd = [0]*trpsNum
    trpsFxd[i] = 1
    fxdTrpsMsk = srv.genFixedTrapsMask(trpsFxd)
    mutChrom = srv.mutateChromosome(initChrom, fxdTrpsMsk)
    resSum = np.sum(np.isclose(initChrom, mutChrom))
    results.extend([resSum == dims])
all(results)


trpsFxd = [0]*trpsNum
results = []
for i in range(trpsNum):
    trpsFxd[i] = 1
    fxdTrpsMsk = srv.genFixedTrapsMask(trpsFxd)
    mutChrom = srv.mutateChromosome(initChrom, fxdTrpsMsk)
    resSum = np.sum(np.isclose(initChrom, mutChrom))
    results.extend([resSum])