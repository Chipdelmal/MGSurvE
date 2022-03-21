
import math
import MGSurvE as srv

# Path for outputs
out_pth = './sims_out/'
# Landscape's bounding box
bbox = ((-350, 350), (-225, 225))
# Mosquito movement kernel
mKer = {'params': [.075, 1.0e-10, math.inf], 'zeroInflation': .75}
# Number of sites and clusters in the environment
ptsNum = 350
(clsNum, clsRad) = (5, 75)
# Probability for each point-type
pTypesProb =[0.05, 0.70, 0.25]
# Number and type of traps
nullTraps = [0, 0, 0, 0]
typeTraps = [0, 0, 1, 1]
# Traps' kernels
tKer = {
    0: {'kernel': srv.exponentialDecay, 'params': {'A': .5, 'b': .075}},
    1: {'kernel': srv.exponentialDecay, 'params': {'A': .25, 'b': .025}}
}
# Transition probabilities between point-types
msk = [
    [0.05, 0.70, 0.25],
    [0.30, 0.10, 0.60],
    [0.70, 0.10, 0.20],
]
# GA Settings
(gens, verbose) = (2500, False)
gaParams = [
    {'mate': .3, 'cxpb': 0.5}, 
    {'mean': 0, 'sd': min([i[1]-i[0] for i in bbox])/5, 'mutpb': .4, 'ipb': .5},
    {'tSize': 3}
]
# Plots
(dpi, pad) = (350, 0.05)
