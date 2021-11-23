
import math
import unittest
import numpy as np
import pandas as pd
import MGSurvE as srv

###############################################################################
# Define landscape for tests
###############################################################################
points = pd.DataFrame({
    'x': [0, 0, 1],
    'y': [0, 2, 1],
    't': [0, 1, 0]
})
msk = [
    [.9, .1],
    [.1, .9]
]
# Traps info ------------------------------------------------------------------
traps = pd.DataFrame({
    'x': [0, 1],
    'y': [2, 2],
    't': ['b', 'a']
})
tker = {
    'a': {'kernel': srv.exponentialDecay, 'params': srv.BASIC_EXP_TRAP},
    'b': {'kernel': srv.exponentialDecay, 'params': {'A': 0.1, 'b': 0.5}} 
}
# Generating landscape --------------------------------------------------------
lnd = srv.Landscape(
    points, maskingMatrix=msk, 
    traps=traps, trapsKernels=tker
)


def test_SelectiveMutation():
    # Fix all traps -----------------------------------------------------------
    trapsNew = pd.DataFrame({
        'x': [0, 0, 1, 1], 'y': [1, 1, 0, 0],
        't': [1, 0, 1, 0], 'f': [1, 1, 1, 1]
    })
    tkerNew = {
        0: {'kernel': srv.exponentialDecay, 'params': {'A': 100000, 'b': 0}},
        1: {'kernel': srv.exponentialDecay, 'params': {'A': 0, 'b': 0}} 
    }
    lnd.updateTraps(trapsNew, tkerNew)
    chrom = srv.initChromosome(lnd.trapsNumber, coordsRange=((0, 5), (0, 5)))
    fxdTrpsMsk = srv.genFixedTrapsMask(lnd.trapsFixed)
    mut = srv.mutateChromosome(chrom, fxdTrpsMsk)
    noMutation = all(np.equal(chrom, mut))
    # Move all traps ----------------------------------------------------------
    trapsNew = pd.DataFrame({
        'x': [0, 0, 1, 1], 'y': [1, 1, 0, 0],
        't': [1, 0, 1, 0], 'f': [0, 0, 0, 0]
    })
    tkerNew = {
        0: {'kernel': srv.exponentialDecay, 'params': {'A': 100000, 'b': 0}},
        1: {'kernel': srv.exponentialDecay, 'params': {'A': 0, 'b': 0}} 
    }
    lnd.updateTraps(trapsNew, tkerNew)
    chrom = srv.initChromosome(lnd.trapsNumber, coordsRange=((0, 5), (0, 5)))
    fxdTrpsMsk = srv.genFixedTrapsMask(lnd.trapsFixed)
    mut = srv.mutateChromosome(chrom, fxdTrpsMsk)
    mutation = any(np.equal(chrom, mut))
    # Move one trap -----------------------------------------------------------
    trapsNew = pd.DataFrame({
        'x': [0, 0, 1, 1], 'y': [1, 1, 0, 0],
        't': [1, 0, 1, 0], 'f': [0, 1, 0, 0]
    })
    tkerNew = {
        0: {'kernel': srv.exponentialDecay, 'params': {'A': 100000, 'b': 0}},
        1: {'kernel': srv.exponentialDecay, 'params': {'A': 0, 'b': 0}} 
    }
    lnd.updateTraps(trapsNew, tkerNew)
    chrom = srv.initChromosome(lnd.trapsNumber, coordsRange=((0, 5), (0, 5)))
    fxdTrpsMsk = srv.genFixedTrapsMask(lnd.trapsFixed)
    mut = srv.mutateChromosome(chrom, fxdTrpsMsk)
    mutationOne = any(np.equal(chrom, mut))
    # Combine tests -----------------------------------------------------------
    assert (noMutation and not(mutation) and mutationOne)


def test_selectiveCrossover():
    # Fix all traps -----------------------------------------------------------
    trapsNew = pd.DataFrame({
        'x': [0, 0, 1, 1], 'y': [1, 1, 0, 0],
        't': [0, 0, 0, 0], 'f': [1, 1, 1, 1]
    })
    tkerNew = {
        0: {'kernel': srv.exponentialDecay, 'params': {'A': 1, 'b': 0.5}}
    }
    lnd.updateTraps(trapsNew, tkerNew)
    chrom = srv.initChromosome(lnd.trapsNumber, coordsRange=((0, 5), (0, 5)))
    fxdTrpsMsk = srv.genFixedTrapsMask(lnd.trapsFixed)
    mut = srv.mutateChromosome(chrom, fxdTrpsMsk)
    noMutation = all(np.equal(chrom, mut))
    (pre1, pre2) = (chrom.copy(), mut.copy())
    (ind1, ind2) = srv.cxBlend(chrom, mut, fxdTrpsMsk)
    fxdA = all([np.isclose(a, b) for (a, b) in zip(pre1, ind1)])
    fxdB = all([np.isclose(a, b) for (a, b) in zip(pre2, ind2)])
    # Move all traps ----------------------------------------------------------
    trapsNew = pd.DataFrame({
        'x': [0, 0, 1, 1], 'y': [1, 1, 0, 0],
        't': [0, 0, 0, 0], 'f': [0, 0, 0, 0]
    })
    tkerNew = {
        0: {'kernel': srv.exponentialDecay, 'params': {'A': 1, 'b': 0.5}}
    }
    lnd.updateTraps(trapsNew, tkerNew)
    chrom = srv.initChromosome(lnd.trapsNumber, coordsRange=((0, 5), (0, 5)))
    fxdTrpsMsk = srv.genFixedTrapsMask(lnd.trapsFixed)
    mut = srv.mutateChromosome(chrom, fxdTrpsMsk)
    (pre1, pre2) = (chrom.copy(), mut.copy())
    (ind1, ind2) = srv.cxBlend(chrom, mut, fxdTrpsMsk)
    movA = any([np.isclose(a, b) for (a, b) in zip(pre1, ind1)])
    movB = any([np.isclose(a, b) for (a, b) in zip(pre2, ind2)])
    # Combine tests -----------------------------------------------------------
    assert (fxdA and fxdB) and (not (movA or movB))



###############################################################################
# Main
###############################################################################
if __name__ == '__main__':
    unittest.main()