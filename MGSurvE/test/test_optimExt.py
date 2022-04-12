
import math
import unittest
import numpy as np
import pandas as pd
from copy import deepcopy
import MGSurvE as srv


def test_initChromosomeMixedDataType():
    # Setting landscape up ----------------------------------------------------
    pts = ((-100, -50, 0), (100, 50, 0))
    points = pd.DataFrame(pts, columns=('x', 'y', 't'))
    (nullC, ktyp) = ([0, 0, 0, 0, 0, 0], [1, 3, 2, 1, 0, 1])
    traps = pd.DataFrame({
        'x': nullC, 'y': nullC, 't': ktyp,
        'f': [1, 1, 1, 1, 1, 1],
        'o': [0, 0, 0, 0, 0, 0]
    })
    tKer = {
        0: {'kernel': srv.exponentialDecay, 'params': {'A': .3, 'b': .05}},
        1: {'kernel': srv.exponentialDecay, 'params': {'A': .35, 'b': .04}},
        2: {'kernel': srv.exponentialDecay, 'params': {'A': .25,  'b': .025}} ,
        3: {'kernel': srv.sigmoidDecay,     'params': {'A': .2, 'rate': .5, 'x0': 1}}
    }
    # Generate landscape ------------------------------------------------------
    lnd = srv.Landscape(points, traps=traps, trapsKernels=tKer)
    bbox = lnd.getBoundingBox()
    trpMsk = srv.genFixedTrapsMask(lnd.trapsFixed)
    trpTsk = lnd.trapsTOptim
    # Init chromosome ---------------------------------------------------------
    chromBase = srv.initChromosomeMixed(
        trapsCoords=lnd.trapsCoords, 
        fixedTrapsMask=trpMsk, typeOptimMask=trpTsk,
        coordsRange=bbox, indpb=1,
        trapsPool=[0, 1, 2, 3, 0, 1, 2]   
    )
    # Test for type -----------------------------------------------------------
    results = [False]*len(chromBase)
    for (ix, val) in enumerate(chromBase):
        if ix < (traps.shape[0]*2):
            results[ix] = isinstance(val, float)
        else:
            results[ix] = isinstance(val, int)
    # Combine tests -----------------------------------------------------------
    assert all(results)


def test_initChromosomeMixedCoords():
    # Setting landscape up ----------------------------------------------------
    pts = ((-100, -50, 0), (100, 50, 0))
    points = pd.DataFrame(pts, columns=('x', 'y', 't'))
    (nullC, ktyp) = ([0, 0, 0, 0, 0, 0], [1, 3, 2, 1, 0, 1])
    traps = pd.DataFrame({
        'x': nullC, 'y': nullC, 't': ktyp,
        'f': [1, 1, 1, 1, 1, 1],
        'o': [0, 0, 0, 0, 0, 0]
    })
    tKer = {
        0: {'kernel': srv.exponentialDecay, 'params': {'A': .3, 'b': .05}},
        1: {'kernel': srv.exponentialDecay, 'params': {'A': .35, 'b': .04}},
        2: {'kernel': srv.exponentialDecay, 'params': {'A': .25,  'b': .025}} ,
        3: {'kernel': srv.sigmoidDecay,     'params': {'A': .2, 'rate': .5, 'x0': 1}}
    }
    lnd = srv.Landscape(points, traps=traps, trapsKernels=tKer)
    lndBase = deepcopy(lnd)
    bbox = lnd.getBoundingBox()
    # Init Baseline -----------------------------------------------------------
    trapsPool = list(traps['t'])+list(range(100, 300))
    trpsNum = traps.shape[0]
    baseChrom = srv.initChromosomeMixed(
        trapsCoords=lndBase.trapsCoords, 
        fixedTrapsMask=srv.genFixedTrapsMask(lndBase.trapsFixed), 
        typeOptimMask=lnd.trapsTOptim,
        coordsRange=bbox, indpb=1,
        trapsPool=trapsPool
    )
    # Hand-test coords mutation -----------------------------------------------
    lndTest = deepcopy(lndBase)
    traps = pd.DataFrame({
        'x': nullC, 'y': nullC, 't': ktyp,
        'f': [0, 1, 0, 1, 0, 1],
        'o': [0, 0, 0, 0, 0, 0]
    })
    lndTest.updateTraps(traps, tKer)
    testChrom = srv.initChromosomeMixed(
        trapsCoords=lndTest.trapsCoords, 
        fixedTrapsMask=srv.genFixedTrapsMask(lndTest.trapsFixed), 
        typeOptimMask=lndTest.trapsTOptim,
        coordsRange=bbox, indpb=1,
        trapsPool=trapsPool
    )
    # Check if half the coordinates did mutate
    coordsSect = [i[:trpsNum*2] for i in (baseChrom, testChrom)]
    passed = sum([np.isclose(a, b) for (a, b) in zip(*coordsSect)])
    coordsPass = (passed == (trpsNum*2)/2)
    # Check if types optim stayed consistent
    typesSect = [i[trpsNum*2:] for i in (baseChrom, testChrom)]
    passed = [a == b for (a, b) in zip(*typesSect)]
    typesPass = all(passed)
    # Put tests together ------------------------------------------------------
    assert (coordsPass and typesPass)


def test_initChromosomeMixedTypes():
    # Setting landscape up ----------------------------------------------------
    pts = ((-100, -50, 0), (100, 50, 0))
    points = pd.DataFrame(pts, columns=('x', 'y', 't'))
    (nullC, ktyp) = ([0, 0, 0, 0, 0, 0], [1, 3, 2, 1, 0, 1])
    traps = pd.DataFrame({
        'x': nullC, 'y': nullC, 't': ktyp,
        'f': [1, 1, 1, 1, 1, 1],
        'o': [0, 0, 0, 0, 0, 0]
    })
    tKer = {
        0: {'kernel': srv.exponentialDecay, 'params': {'A': .3, 'b': .05}},
        1: {'kernel': srv.exponentialDecay, 'params': {'A': .35, 'b': .04}},
        2: {'kernel': srv.exponentialDecay, 'params': {'A': .25,  'b': .025}} ,
        3: {'kernel': srv.sigmoidDecay,     'params': {'A': .2, 'rate': .5, 'x0': 1}}
    }
    lnd = srv.Landscape(points, traps=traps, trapsKernels=tKer)
    lndBase = deepcopy(lnd)
    bbox = lnd.getBoundingBox()
    # Init Baseline -----------------------------------------------------------
    trapsPool = list(traps['t'])+list(range(100, 300))
    trpsNum = traps.shape[0]
    baseChrom = srv.initChromosomeMixed(
        trapsCoords=lndBase.trapsCoords, 
        fixedTrapsMask=srv.genFixedTrapsMask(lndBase.trapsFixed), 
        typeOptimMask=lnd.trapsTOptim,
        coordsRange=bbox, indpb=1,
        trapsPool=trapsPool
    )
    # Hand-test coords mutation -----------------------------------------------
    lndTest = deepcopy(lndBase)
    traps = pd.DataFrame({
        'x': nullC, 'y': nullC, 't': ktyp,
        'f': [1, 1, 1, 1, 1, 1],
        'o': [0, 1, 0, 1, 0, 1]
    })
    lndTest.updateTraps(traps, tKer)
    testChrom = srv.initChromosomeMixed(
        trapsCoords=lndTest.trapsCoords, 
        fixedTrapsMask=srv.genFixedTrapsMask(lndTest.trapsFixed), 
        typeOptimMask=lndTest.trapsTOptim,
        coordsRange=bbox, indpb=1,
        trapsPool=trapsPool
    )
    # Check if half the coordinates did not mutate
    coordsSect = [i[:trpsNum*2] for i in (baseChrom, testChrom)]
    passed = sum([np.isclose(a, b) for (a, b) in zip(*coordsSect)])
    coordsPass = (passed == (trpsNum*2))
    # Check if types section was mutated 
    typesSect = [i[trpsNum*2:trpsNum*2+trpsNum] for i in (baseChrom, testChrom)]
    passed = sum([a == b for (a, b) in zip(*typesSect)])
    typesPass = (passed > trpsNum*.4)
    print(typesPass)
    # Put tests together ------------------------------------------------------
    assert (coordsPass and typesPass)

###############################################################################
# Main
###############################################################################
if __name__ == '__main__':
    unittest.main()