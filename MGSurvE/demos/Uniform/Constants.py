
import math
import MGSurvE as srv

# Path for outputs
out_pth = './sims_out/'
# Landscape's bounding box
bbox = ((-100, 100), (-100, 100))
# Mosquito movement kernel
mKerZ = mKer = {'params': srv.MEDIUM_MOV_EXP_PARAMS, 'zeroInflation': .50} 
mKerN = mKer = {'params': srv.MEDIUM_MOV_EXP_PARAMS, 'zeroInflation': .25} 
# Number of sites and clusters in the environment
ptsNum = 150
(clsNum, clsRad) = (5, 175)
# Probability for each point-type
pTypesProb =[0.1, 0.7, 0.2]
# Number and type of traps
nullTraps = [0, 0, 0, 0, 0]
typeTraps = [0, 0, 0, 1, 0]
# Traps' kernels
tKer = {
    0: {'kernel': srv.exponentialDecay, 'params': {'A': .75, 'b': .100}},
    1: {'kernel': srv.exponentialDecay, 'params': {'A': .50, 'b': .125}}
}
# Transition probabilities between point-types
msk = [
    [0.05, 0.90, 0.05],
    [0.25, 0.05, 0.70],
    [0.90, 0.05, 0.05],
]
# GA Settings
(gens, verbose) = (10000, False)
gaParams = [
    {
        'mate': .3, 
        'cxpb': 0.3,
        'indpb': 0.5,
        'alpha': .5
    }, 
    {
        'mean': 0, 'sd': min([i[1]-i[0] for i in bbox])/5, 
        'mutpb': .4, 'ipb': .5, 
        'indpb': 0.5
    },
    {'tSize': 4}
]
# Plots
(dpi, pad, pad_i) = (350, 0.05, (10, 10))
