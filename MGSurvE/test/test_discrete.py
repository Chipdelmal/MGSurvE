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


###############################################################################
# Main
###############################################################################
if __name__ == '__main__':
    unittest.main()