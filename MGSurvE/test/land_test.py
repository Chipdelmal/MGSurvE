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
    points, maskingMatrix=msk, traps=traps, trapsKernels=tker
)

###############################################################################
# Tests on landscape
###############################################################################
def test_LandPointsNumbers():
    ptsMatch = (lnd.pointNumber == points.shape[0])
    trpMatch = (lnd.trapsNumber == traps.shape[0])
    full = all([ptsMatch, trpMatch])
    assert full

def test_MarkovMatrix():


if __name__ == '__main__':
    unittest.main()