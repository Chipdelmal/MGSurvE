
import unittest
import numpy as np
from copy import deepcopy
import MGSurvE as srv

def test_InitDiscrete():
    (vct, ban, ids) = ([0]*100, {5, 6}, 25)
    chrom = srv.initDiscreteChromosome(range(ids), vct, ban)
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
    trpsNum = 25
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

###############################################################################
# Main
###############################################################################
if __name__ == '__main__':
    unittest.main()