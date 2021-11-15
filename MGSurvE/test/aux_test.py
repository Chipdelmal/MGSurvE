
import math
import unittest
import numpy as np
import pandas as pd
import MGSurvE as srv


###############################################################################
# Tests on landscape
###############################################################################
def test_kSolve():
    kDict = {'kernel': srv.exponentialDecay, 'params': srv.BASIC_EXP_TRAP}
    sols = [
        srv.nSolveKernel({
            'kernel': srv.exponentialDecay, 
            'params': {'A': 0.5, 'b': 0.15}
        }, i)
        for i in [0.5, 0.00000001]
    ]
    assert (np.isclose(sols[0], 0) and (sols[1] > 100))

###############################################################################
# Main
###############################################################################
if __name__ == '__main__':
    unittest.main()