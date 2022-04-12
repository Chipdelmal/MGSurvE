
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
    traps = pd.DataFrame({
        'x': [0, 0, 0, 0, 0, 0],
        'y': [0, 0, 0, 0, 0, 0],
        't': [1, 3, 2, 1, 0, 1],
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


###############################################################################
# Main
###############################################################################
if __name__ == '__main__':
    unittest.main()