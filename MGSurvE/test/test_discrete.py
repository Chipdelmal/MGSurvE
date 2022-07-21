import unittest
import MGSurvE as srv

def test_InitDiscrete():
    (vct, ban, ids) = ([0]*100, {5, 6}, 25)
    chrom = srv.initDiscreteChromosome(range(ids), vct, ban)
    vix = set(chrom)
    # Check for conditions ----------------------------------------------------
    lTest = (len(chrom) == len(vct))
    ixTest = ((vix-ban)==vix)
    assert (lTest and ixTest)

###############################################################################
# Main
###############################################################################
if __name__ == '__main__':
    unittest.main()