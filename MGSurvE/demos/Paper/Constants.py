
import math
import MGSurvE as srv

# Path for outputs
out_pth = './sims_out/'
# Landscape's bounding box
bbox = ((-150, 150), (-150, 150))
# Mosquito movement kernel
mKerZ = {'params': [0.01848777, 1.0e-10, math.inf], 'zeroInflation': .75}
mKerN = {'params': [0.01848777, 1.0e-10, math.inf], 'zeroInflation': 0}
# Number of sites and clusters in the environment
ptsNum = 250
(clsNum, clsRad) = (5, 175)
# Probability for each point-type
pTypesProb =[0.1, 0.7, 0.2]
# Number and type of traps
nullTraps = [0, 0, 0, 0, 0, 0, 0, 0]
typeTraps = [0, 0, 0, 0, 1, 1, 1, 1]
# Traps' kernels
# tKer = {
#     0: {'kernel': srv.exponentialDecay, 'params': {'A': 0.80, 'b': .035}},
#     1: {'kernel': srv.exponentialDecay, 'params': {'A': 1.00, 'b': .040}},
#     2: {'kernel': srv.exponentialDecay, 'params': {'A': 0.75, 'b': .025}}
# }
# Long ------------------------------------------------------------------------
tKer = {
    1: {
        'kernel': srv.sigmoidDecay,     
        'params': {'A': 1.0, 'rate': .25, 'x0': 1/0.0409522}
    },
    0: {
        'kernel': srv.exponentialDecay, 
        'params': {'A': 1.0, 'b': 0.0409522}
    }
}
# Short -----------------------------------------------------------------------
# tKer = {
#     1: {
#         'kernel': srv.sigmoidDecay,     
#         'params': {'A': 1, 'rate': .25, 'x0': 1/0.0629534}
#     },
#     0: {
#         'kernel': srv.exponentialDecay, 
#         'params': {'A': 1, 'b': 0.0629534}
#     }
# }
# Transition probabilities between point-types
msk = [
    [0.05, 0.90, 0.05],
    [0.25, 0.05, 0.70],
    [0.90, 0.05, 0.05],
]
# GA Settings
(gens, verbose) = (1000, True)
gaParams = [
    {
        'mate': .35, 
        'cxpb': 0.5,
        'indpb': 0.35,
        'alpha': .5
    }, 
    {
        'mean': 0, 'sd': min([i[1]-i[0] for i in bbox])/5, 
        'mutpb': .5, 'ipb': .5, 
        'indpb': 0.35
    },
    {'tSize': 5}
]
# Plots
(dpi, pad, pad_i) = (350, 0.05, (10, 10))
