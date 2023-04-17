
import math
import MGSurvE as srv

# Path for outputs
out_pth = './sims_out/'
# Landscape's bounding box
bbox = ((-75, 75), (-75, 75))
# Mosquito movement kernel
mKerZ = {'params': [0.01848777, 1.0e-10, math.inf], 'zeroInflation': .25}
mKerN = {'params': [0.01848777, 1.0e-10, math.inf], 'zeroInflation': 0}
# Number of sites and clusters in the environment
ptsNum = 225
(clsNum, clsRad) = (5, 175)
# Probability for each point-type
pTypesProb =[0.1, 0.7, 0.2]
# Number and type of traps
nullTraps = [0, 0, 0, 0, 0]
typeTraps = [0, 0, 0, 0, 1]
# Traps' kernels
tKer = {
    1: {
        'kernel': srv.sigmoidDecay,     
        'params': {'A': 1, 'rate': .25, 'x0': 1/0.02}
    },
    0: {
        'kernel': srv.exponentialDecay, 
        'params': {'A': 1, 'b': 0.02}
    }
}
# Transition probabilities between point-types
msk = [
    [0.05, 0.90, 0.05],
    [0.25, 0.05, 0.70],
    [0.90, 0.05, 0.05],
]
# GA Settings
(gens, verbose) = (5000, False)
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
