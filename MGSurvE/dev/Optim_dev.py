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
# Crossover function test
###############################################################################


(trpsNum, dims) = (12, 2)
trpsFxd = [0]*trpsNum
result = []
for i in range(trpsNum):
    trpsFxd = [0]*trpsNum
    trpsFxd[i] = 1
    fxdTrpsMsk = srv.genFixedTrapsMask(trpsFxd)
    chromA = srv.initChromosome(trpsNum, coordsRange=((0, 5), (0, 5)))
    chromB = srv.initChromosome(trpsNum, coordsRange=((0, 5), (0, 5)))
    (pre1, pre2) = (chromA.copy(), chromB.copy())
    (ind1, ind2) = srv.cxBlend(chromA, chromB, fxdTrpsMsk)
    fxdA = (np.sum([np.isclose(a, b) for (a, b) in zip(pre1, ind1)]) == dims)
    fxdB = (np.sum([np.isclose(a, b) for (a, b) in zip(pre2, ind2)]) == dims)
    result.extend([fxdA and fxdB])
testShift = all(result)



trpsFxd = [0]*trpsNum
(total, result) = (0, [])
for i in range(trpsNum):
    total = total + dims
    trpsFxd[i] = 1
    fxdTrpsMsk = srv.genFixedTrapsMask(trpsFxd)
    chromA = srv.initChromosome(trpsNum, coordsRange=((0, 5), (0, 5)))
    chromB = srv.initChromosome(trpsNum, coordsRange=((0, 5), (0, 5)))
    (pre1, pre2) = (chromA.copy(), chromB.copy())
    (ind1, ind2) = srv.cxBlend(chromA, chromB, fxdTrpsMsk)
    fxdA = np.sum([np.isclose(a, b) for (a, b) in zip(pre1, ind1)])
    fxdB = np.sum([np.isclose(a, b) for (a, b) in zip(pre2, ind2)])
    result.extend([fxdA == fxdB == total])
testCumsum = all(result)