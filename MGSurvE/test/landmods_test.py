
import math
import unittest
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import MGSurvE as srv


###############################################################################
# Generate landscape object
###############################################################################
ptsNum = 20
bbox = ((-10, 10), (-10, 10))
xy = srv.ptsRandUniform(ptsNum, bbox).T
points = pd.DataFrame({'x': xy[0], 'y': xy[1], 't': [0]*xy.shape[1]})
lndOrg = srv.Landscape(points, attractionVector=[1]*ptsNum)
orgMig = lndOrg.maskedMigration.T

def test_LandAttractiveness():
    (HI, LO) = (100, 5)
    lndMod = srv.Landscape(
        points, attractionVector=([HI]*(ptsNum//2)) +([LO]*(ptsNum//2))
    )
    modMig = lndMod.maskedMigration.T
    # Test Higher-Level -------------------------------------------------------
    difference = not(all(np.isclose(np.sum(orgMig, axis=1), np.sum(modMig, axis=1))))
    consistency = all(np.isclose(np.sum(orgMig, axis=0), np.sum(modMig, axis=0)))
    # Test I/O ----------------------------------------------------------------
    (diff, ixSplt) = ((modMig.T - orgMig.T), ptsNum//2)
    # Hi In-degree (all iHi are higher than original)
    hiIn = sum([sum(i) for i in (diff[:,:ixSplt] > 0)])
    # Lo In-degree (all oLo are lower than original)
    loIn = sum([sum(i) for i in (diff[:,ixSplt:] < 0)])
    threshold = [i > (ptsNum*ixSplt*.9) for i in (hiIn, loIn)]
    # Assert tests 
    assert all([consistency, difference]+threshold)



###############################################################################
# Main
###############################################################################
if __name__ == '__main__':
    unittest.main()