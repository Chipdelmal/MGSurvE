
import math
import unittest
import numpy as np
import pandas as pd
import MGSurvE as srv


# def test_SelectiveMutation():
#     (trpsNum, dims) = (12, 2)
#     trpsFxd = [0]*trpsNum
#     initChrom = srv.initChromosome(trpsNum, coordsRange=((0, 5), (0, 5)))
#     # Test one shift mask -----------------------------------------------------
#     results = []
#     for i in range(trpsNum):
#         trpsFxd = [0]*trpsNum
#         trpsFxd[i] = 1
#         fxdTrpsMsk = srv.genFixedTrapsMask(trpsFxd)
#         mutChrom = srv.mutateChromosome(initChrom, fxdTrpsMsk)
#         resSum = np.sum(np.isclose(initChrom, mutChrom))
#         results.extend([resSum == dims])
#     testShift = all(results)
#     # Test cumulative on mask -------------------------------------------------
#     (trpsFxd, results, total) = ([0]*trpsNum, [], 0)
#     for i in range(trpsNum):
#         total = total + (i+1)
#         trpsFxd[i] = 1
#         fxdTrpsMsk = srv.genFixedTrapsMask(trpsFxd)
#         mutChrom = srv.mutateChromosome(initChrom, fxdTrpsMsk)
#         resSum = np.sum(np.isclose(initChrom, mutChrom))
#         results.extend([resSum])
#     testCumsum = (np.sum(results)//2 == total)
#     # Combine tests -----------------------------------------------------------
#     assert (testShift and testCumsum)


# def test_selectiveCrossover():
#     (trpsNum, dims) = (12, 2)
#     trpsFxd = [0]*trpsNum
#     result = []
#     # Test one shift mask -----------------------------------------------------
#     for i in range(trpsNum):
#         trpsFxd = [0]*trpsNum
#         trpsFxd[i] = 1
#         fxdTrpsMsk = srv.genFixedTrapsMask(trpsFxd)
#         chromA = srv.initChromosome(trpsNum, coordsRange=((0, 5), (0, 5)))
#         chromB = srv.initChromosome(trpsNum, coordsRange=((0, 5), (0, 5)))
#         (pre1, pre2) = (chromA.copy(), chromB.copy())
#         (ind1, ind2) = srv.cxBlend(chromA, chromB, fxdTrpsMsk)
#         fxdA = (np.sum([np.isclose(a, b) for (a, b) in zip(pre1, ind1)]) == dims)
#         fxdB = (np.sum([np.isclose(a, b) for (a, b) in zip(pre2, ind2)]) == dims)
#         result.extend([fxdA and fxdB])
#     testShift = all(result)
#     # Test cumulative on mask -------------------------------------------------
#     trpsFxd = [0]*trpsNum
#     (total, result) = (0, [])
#     for i in range(trpsNum):
#         total = total + dims
#         trpsFxd[i] = 1
#         fxdTrpsMsk = srv.genFixedTrapsMask(trpsFxd)
#         chromA = srv.initChromosome(trpsNum, coordsRange=((0, 5), (0, 5)))
#         chromB = srv.initChromosome(trpsNum, coordsRange=((0, 5), (0, 5)))
#         (pre1, pre2) = (chromA.copy(), chromB.copy())
#         (ind1, ind2) = srv.cxBlend(chromA, chromB, fxdTrpsMsk)
#         fxdA = np.sum([np.isclose(a, b) for (a, b) in zip(pre1, ind1)])
#         fxdB = np.sum([np.isclose(a, b) for (a, b) in zip(pre2, ind2)])
#         result.extend([fxdA == fxdB == total])
#     testCumsum = all(result)
#     # Combine tests -----------------------------------------------------------
#     assert (testShift and testCumsum)


###############################################################################
# Main
###############################################################################
if __name__ == '__main__':
    unittest.main()