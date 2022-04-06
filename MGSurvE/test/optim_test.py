
import math
import unittest
import numpy as np
import pandas as pd
from copy import deepcopy
import MGSurvE as srv

def test_SelectiveMutation_OneShift():
    (trpsNum, dims) = (12, 2)
    trpsFxd = [0]*trpsNum
    initMsk = [True] * trpsNum * dims
    chrom = np.reshape(np.random.uniform(-10, 10, trpsNum*2), (-1, dims))
    initChrom = srv.initChromosome(chrom, initMsk, coordsRange=((-10, 10), (-10, 10)))
    # Test one shift mask -----------------------------------------------------
    results = []
    i = 0
    for i in range(trpsNum):
        trpsFxd = [0]*trpsNum
        trpsFxd[i] = 1
        fxdTrpsMsk = srv.genFixedTrapsMask(trpsFxd)
        mutChrom = deepcopy(initChrom)
        srv.mutateChromosome(mutChrom, fxdTrpsMsk, randArgs={'loc': 10}, indpb=1)
        resSum = np.sum(np.isclose(initChrom, mutChrom))
        results.extend([resSum == dims])
    testShift = all(results)
    assert testShift


def test_SelectiveMutation_CumShift():
    (trpsNum, dims) = (12, 2)
    trpsFxd = [0]*trpsNum
    initMsk = [True] * trpsNum * dims
    chrom = np.reshape(np.random.uniform(-10, 10, trpsNum*2), (-1, dims))
    initChrom = srv.initChromosome(chrom, initMsk, coordsRange=((-10, 10), (-10, 10)))
    # Test cumulative on mask -------------------------------------------------
    (trpsFxd, results, total) = ([0]*trpsNum, [], 0)
    for i in range(trpsNum):
        total = total + (i+1)
        trpsFxd[i] = 1
        fxdTrpsMsk = srv.genFixedTrapsMask(trpsFxd)
        mutChrom = deepcopy(initChrom)
        srv.mutateChromosome(initChrom, fxdTrpsMsk, randArgs={'loc': 10}, indpb=1)
        resSum = np.sum(np.isclose(initChrom, mutChrom))
        results.extend([resSum])
    testCumsum = (np.sum(results)//2 == total)
    assert testCumsum


def test_selectiveCrossover_OneShift():
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
    assert testShift


def test_selectiveCrossover_CumShift():
    (trpsNum, dims) = (12, 2)
    trpsFxd = [0]*trpsNum
    initMsk = [True] * trpsNum * dims
    result = []
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
    assert testCumsum


def test_mutateChromosomeAsymmetric():
    (trpsNum, dims) = (12, 2)
    trpsFxd = [0]*trpsNum
    initMsk = [True] * trpsNum * dims
    chrom = np.random.uniform(-10, 10, trpsNum*2)
    trpMsk = srv.genFixedTrapsMask(trpsFxd)

    org = np.copy(chrom)
    # Test mutation over X ----------------------------------------------------
    mod = srv.mutateChromosomeAsymmetric(
        np.copy(chrom), trpMsk, 
        randArgs={
            'x': {'loc': 10, 'scale': 100}, 
            'y': {'loc': 0, 'scale': 0}
        },
        indpb=1
    )
    totalX = np.sum([org[a]!=mod[0][a] for a in range(len(org))])
    # Test mutation over Y ----------------------------------------------------
    mod = srv.mutateChromosomeAsymmetric(
        np.copy(chrom), trpMsk, 
        randArgs={
            'x': {'loc': 0, 'scale': 0}, 
            'y': {'loc': 10, 'scale': 100}
        },
        indpb=1
    )
    totalY = np.sum([org[a]!=mod[0][a] for a in range(len(org))])
    # Test mutation over XY --------------------------------------------------
    mod = srv.mutateChromosomeAsymmetric(
        np.copy(chrom), trpMsk, 
        randArgs={
            'x': {'loc': 10, 'scale': 0}, 
            'y': {'loc': 10, 'scale': 1}
        },
        indpb=1
    )
    totalXY = np.sum([org[a]!=mod[0][a] for a in range(len(org))])
    # Combine tests -----------------------------------------------------------
    assert ((totalX==trpsNum) and (totalY==trpsNum) and (totalXY==2*trpsNum))

###############################################################################
# Main
###############################################################################
if __name__ == '__main__':
    unittest.main()