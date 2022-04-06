
import math
import unittest
import numpy as np
import pandas as pd
import MGSurvE as srv
import warnings
warnings.filterwarnings('ignore', 'The iteration is not making good progress')


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
    't': [1, 0]
})
tker = {
    0: {'kernel': srv.exponentialDecay, 'params': srv.BASIC_EXP_TRAP},
    1: {'kernel': srv.exponentialDecay, 'params': {'A': 0.1, 'b': 0.5}} 
}
# Generating landscape --------------------------------------------------------
lnd = srv.Landscape(
    points, maskingMatrix=msk, 
    traps=traps, trapsKernels=tker
)

###############################################################################
# Tests on landscape
###############################################################################
def test_LandPointsNumbers():
    ptsMatch = (lnd.pointNumber == points.shape[0])
    trpMatch = (lnd.trapsNumber == traps.shape[0])
    full = all([ptsMatch, trpMatch])
    assert full

def test_MarkovMatrices():
    tsts = []
    for mat in (lnd.migrationMatrix, lnd.maskedMigration, lnd.trapsMigration):
        sumsOne = all([np.isclose(i, 1) for i in np.sum(mat, axis=1)])
        tsts.extend([sumsOne])
    assert all(tsts)

def test_UpdateMigration():
    trapsNew = pd.DataFrame({
        'x': [0, 0],
        'y': [1, 1],
        't': [1, 0],
        'f': [0, 0]
    })
    tkerNew = {
        0: {'kernel': srv.exponentialDecay, 'params': {'A': 100000, 'b': 0}},
        1: {'kernel': srv.exponentialDecay, 'params': {'A': 0, 'b': 0}} 
    }
    lnd.updateTraps(trapsNew, tkerNew)
    sumsNumber = lnd.pointNumber+lnd.trapsNumber-1
    maxTrap = np.isclose(np.sum(lnd.trapsMigration[:,-1]), sumsNumber)
    nullTrap = np.isclose(np.sum(lnd.trapsMigration[:,-2]), 1)
    test_MarkovMatrices()
    assert(all([maxTrap, nullTrap]))

def test_MarkovFundamentalMatrix():
    (tau, sitesN, trapsN) = (
        lnd.trapsMigration, lnd.pointNumber, lnd.trapsNumber
    )
    tauC = srv.reshapeInCanonicalForm(tau, lnd.pointNumber, lnd.trapsNumber)
    F_A = srv.getMarkovAbsorbing(tauC, trapsN)
    F_B = srv.getFundamentalMatrix(tau, sitesN, trapsN)
    assert np.equal(F_A, F_B).all()


def test_TrapTypeIndexMask():
    traps = pd.DataFrame({
        'x': [0, 1, 2, 0],
        'y': [2, 2, 0, 1],
        't': [1, 1, 1, 2]
    })
    tker = {
        0: {'kernel': srv.exponentialDecay, 'params': srv.BASIC_EXP_TRAP},
        1: {'kernel': srv.exponentialDecay, 'params': {'A': 0.1, 'b': 0.5}},
        2: {'kernel': srv.exponentialDecay, 'params': {'A': 0.1, 'b': 0.5}} 
    }
    # Generating landscape --------------------------------------------------------
    lnd = srv.Landscape(
        points, maskingMatrix=msk, 
        traps=traps, trapsKernels=tker
    )
    shapes = (lnd.trapsMask.shape == (max(traps['t'])+1, len(set(points['t']))))
    blank = all([all(i==1) for i in lnd.trapsMask])
    assert all([shapes, blank])

###############################################################################
# Main
###############################################################################
if __name__ == '__main__':
    unittest.main()



