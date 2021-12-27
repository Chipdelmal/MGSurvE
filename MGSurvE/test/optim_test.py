
import math
import unittest
import numpy as np
import pandas as pd
from copy import deepcopy
import MGSurvE as srv


def test_SelectiveMutation():
    (trpsNum, dims) = (12, 2)
    trpsFxd = [0]*trpsNum
    initMsk = [True] * trpsNum * dims
    chrom = np.reshape(np.random.uniform(-10, 10, trpsNum*2), (-1, dims))
    initChrom = srv.initChromosome(chrom, initMsk, coordsRange=((-10, 10), (-10, 10)))
    # Test one shift mask -----------------------------------------------------
    results = []
    for i in range(trpsNum):
        trpsFxd = [0]*trpsNum
        trpsFxd[i] = 1
        fxdTrpsMsk = srv.genFixedTrapsMask(trpsFxd)
        mutChrom = deepcopy(initChrom)
        srv.mutateChromosome(mutChrom, fxdTrpsMsk, randArgs={'loc': 10})
        resSum = np.sum(np.isclose(initChrom, mutChrom))
        results.extend([resSum == dims])
    testShift = all(results)
    # Test cumulative on mask -------------------------------------------------
    (trpsFxd, results, total) = ([0]*trpsNum, [], 0)
    for i in range(trpsNum):
        total = total + (i+1)
        trpsFxd[i] = 1
        fxdTrpsMsk = srv.genFixedTrapsMask(trpsFxd)
        mutChrom = deepcopy(initChrom)
        srv.mutateChromosome(initChrom, fxdTrpsMsk, randArgs={'loc': 10})
        resSum = np.sum(np.isclose(initChrom, mutChrom))
        results.extend([resSum])
    testCumsum = (np.sum(results)//2 == total)
    # Combine tests -----------------------------------------------------------
    assert (testShift and testCumsum)


def test_selectiveCrossover():
    (trpsNum, dims) = (12, 2)
    trpsFxd = [0]*trpsNum
    initMsk = [True] * trpsNum * dims
    result = []
    # Test one shift mask -----------------------------------------------------
    for i in range(trpsNum):
        trpsFxd = [0]*trpsNum
        trpsFxd[i] = 1
        fxdTrpsMsk = srv.genFixedTrapsMask(trpsFxd)
        chA = np.reshape(np.random.uniform(-10, 10, trpsNum*2), (-1, dims))
        chB = np.reshape(np.random.uniform(-10, 10, trpsNum*2), (-1, dims))
        chromA = srv.initChromosome(chA, initMsk, coordsRange=((-10, 10), (-10, 10)))
        chromB = srv.initChromosome(chB, initMsk, coordsRange=((-10, 10), (-10, 10)))
        (pre1, pre2) = (chromA.copy(), chromB.copy())
        (ind1, ind2) = srv.cxBlend(chromA, chromB, fxdTrpsMsk)
        fxdA = (np.sum([np.isclose(a, b) for (a, b) in zip(pre1, ind1)]) == dims)
        fxdB = (np.sum([np.isclose(a, b) for (a, b) in zip(pre2, ind2)]) == dims)
        result.extend([fxdA and fxdB])
    testShift = all(result)
    # Test cumulative on mask -------------------------------------------------
    trpsFxd = [0]*trpsNum
    (total, result) = (0, [])
    for i in range(trpsNum):
        total = total + dims
        trpsFxd[i] = 1
        fxdTrpsMsk = srv.genFixedTrapsMask(trpsFxd)
        chA = np.reshape(np.random.uniform(-10, 10, trpsNum*2), (-1, dims))
        chB = np.reshape(np.random.uniform(-10, 10, trpsNum*2), (-1, dims))
        chromA = srv.initChromosome(chA, initMsk, coordsRange=((-10, 10), (-10, 10)))
        chromB = srv.initChromosome(chB, initMsk, coordsRange=((-10, 10), (-10, 10)))
        (pre1, pre2) = (chromA.copy(), chromB.copy())
        (ind1, ind2) = srv.cxBlend(chromA, chromB, fxdTrpsMsk)
        fxdA = np.sum([np.isclose(a, b) for (a, b) in zip(pre1, ind1)])
        fxdB = np.sum([np.isclose(a, b) for (a, b) in zip(pre2, ind2)])
        result.extend([fxdA == fxdB == total])
    testCumsum = all(result)
    # Combine tests -----------------------------------------------------------
    assert (testShift and testCumsum)


###############################################################################
# Main
###############################################################################
if __name__ == '__main__':
    unittest.main()