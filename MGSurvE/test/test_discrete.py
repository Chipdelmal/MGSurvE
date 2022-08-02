
import unittest
import numpy as np
from copy import deepcopy
from math import isclose
import MGSurvE as srv

def test_InitDiscrete():
    (vct, ban, ids) = ([0]*100, {5, 6}, 25)
    chrom = srv.initDiscreteChromosome(
        ptsIds=range(ids), fixedTraps=vct, trapsSiteID=vct, banSites=ban
    )
    vix = set(chrom)
    # Check for conditions ----------------------------------------------------
    lTest = (len(chrom) == len(vct))
    ixTest = ((vix-ban)==vix)
    assert (lTest and ixTest)

def test_MutateDiscrete():
    (ub, chromSize) = (100, 100)
    chromA = srv.mutateDiscreteChromosome(
        [0]*chromSize, range(1, ub), [0]*chromSize, indpb=1
    )[0]
    chromB = srv.mutateDiscreteChromosome(
        [0]*chromSize, range(1, ub), [0]*chromSize, indpb=0
    )[0]
    # Compile tests -----------------------------------------------------------
    noZero = (len([i for i in chromA if i==0]) == 0)
    allZero = (len([i for i in chromB if i==0]) == len(chromB))
    assert (noZero and allZero)

def test_DiscreteSelectiveMutation_OneShift():
    trpsNum = 25
    (trpsFxd, ptdsId) = ([0]*trpsNum, list(range(5000)))
    initChrom = srv.initDiscreteChromosome(ptdsId, trpsFxd, banSites=None)
    # Test one shift mask -----------------------------------------------------
    results = []
    for i in range(trpsNum):
        trpsFxd = [0]*trpsNum
        trpsFxd[i] = 1
        mutChrom = deepcopy(initChrom)
        srv.mutateDiscreteChromosome(
            mutChrom, ptdsId, trpsFxd, 
            indpb=1, banSites=None
        )
        pSum = np.sum([mutChrom[i]==initChrom[i] for i in range(trpsNum)])
        results.append((1 <= pSum) and (pSum < 5))
    assert all(results)

def test_DiscreteSelectiveMutation_CumShift():
    trpsNum = 150
    (trpsFxd, ptdsId) = ([0]*trpsNum, list(range(5000)))
    initChrom = srv.initDiscreteChromosome(ptdsId, trpsFxd, banSites=None)
    # Test cumulative on mask -------------------------------------------------
    results = []
    for i in range(trpsNum):
        trpsFxd[i] = 1
        mutChrom = deepcopy(initChrom)
        srv.mutateDiscreteChromosome(
            mutChrom, ptdsId, trpsFxd, 
            indpb=1, banSites=None
        )
        pSum = np.sum([mutChrom[i]==initChrom[i] for i in range(trpsNum)])
        results.append((i <= pSum) and (pSum < i+5))
    assert all(results)

def test_DiscreteSelectiveCrossover_CumShift():
    trpsNum = 150
    (trpsFxd, ptdsId) = ([0]*trpsNum, list(range(5000)))
    # Test cumulative on mask -------------------------------------------------
    (total, result) = (0, [])
    for i in range(trpsNum):
        total += 1
        trpsFxd[i] = 1
        chromA = srv.initDiscreteChromosome(ptdsId, trpsFxd, ptdsId, banSites=None)
        chromB = srv.initDiscreteChromosome(ptdsId, trpsFxd, ptdsId, banSites=None)
        (pre1, pre2) = (chromA.copy(), chromB.copy())
        (ind1, ind2) = srv.cxDiscreteUniform(chromA, chromB, trpsFxd, indpb=1)
        fxdA = np.sum([isclose(a, b, abs_tol=.25) for (a, b) in zip(pre1, ind1)])
        fxdB = np.sum([isclose(a, b, abs_tol=.25) for (a, b) in zip(pre2, ind2)])
        boolA = isclose(fxdA, total, abs_tol=3)
        boolB = isclose(fxdB, total, abs_tol=3)
        result.append(boolA and boolB)
    assert all(result)

###############################################################################
# Main
###############################################################################
if __name__ == '__main__':
    unittest.main()