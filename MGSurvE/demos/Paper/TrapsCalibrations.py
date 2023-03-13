#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import MGSurvE as srv


distances = [4.572, 9.144, 18.288, 36.576, 54.864, 73.152]
traps = {
    'An. melas': {
        'calves':   [9166, 2712, 1419, 591, 181,  97],
        'co2':      [1245,  505,  195,  83,  54,  26],
        'control':  [ 200,  182,  217, 185, 164, 131]
    },
    'An. other': {
        'calves':   [ 618,  110,   40,   8,   6,   0],
        'co2':      [ 103,   25,   11,   4,   0,   1],
        'control':  [   4,    5,    3,   3,   1,   2]
    },
    'Ae. other': {
        'calves':   [1139,  292,  92,   31,  19,  12],
        'co2':      [ 238,   75,  39,    7,  12,  12],
        'control':  [  10,    9,  16,   14,  10,  18]
    }
}

trapsGeom = {
    'An. melas': {
        'calves':   [(290,4), ( 79,5), ( 41,4), ( 23,5), (  9,6), (  5,3)],
        'co2':      [( 84,9), ( 25,3), ( 11,7), (  5,9), (  4,1), (  1,6)],
        'control':  [( 13,1), ( 10,4), ( 14,3), ( 12,1), ( 11,3), (  8,7)]
    },
    'An. other': {
        'calves':   [( 21,5), (  5,0), (  2,2), (  0,5), (  0,5), (  0  )],
        'co2':      [(  7,1), (  1,6), (  0,7), (  0,5), (  0  ), (  0,1)],
        'control':  [(  0,3), (  0,3), (  0,2), (  0,3), (  0,1), (  0,2)]
    },
    'Ae. other': {
        'calves':   [( 52,3), ( 14,5), (  5,2), (  1,7), (  1,3), (  0,9)],
        'co2':      [( 14,3), (  4,5), (  2,3), (  0,5), (  0,9), (  0,9)],
        'control':  [(  0,7), (  0,7), (  1,5), (  1,1), (  0,8), (  1,4)]
    }
}

spe = 'An. melas'
(fig, ax) = plt.subplots(figsize=(4, 10))
for typ in ['calves', 'co2']:
    sym = ('x' if typ=='calves' else '1')
    trps = [i[0] for i in np.array(trapsGeom[spe][typ])]
    sums = np.sum(trps)
    plt.plot(
        distances, trps/sums*100,
        sym, label=typ
    )
ax.set_yscale('symlog')
ax.set_xlim(0, 80)
ax.set_ylim(0, 100)
plt.legend()

###############################################################################
# Fit Data
###############################################################################
def exponentialDecay(x, B):
    prob = 1 * np.exp(-B * x)
    return prob

def sigmoidDecay(x, rate, x0):
    prob = 1 - 1 / (1 + math.e ** (-rate * (x - x0)))
    return prob

spe = 'An. melas'
(fig, ax) = plt.subplots(figsize=(10, 5))
trps = np.array([i[0] for i in np.array(trapsGeom[spe]['calves'])])+np.array([i[0] for i in np.array(trapsGeom[spe]['co2'])])
sums = np.sum(trps)
plt.plot([0]+distances, [1]+list(trps/sums))
ax.set_xlim(0, 80)
ax.set_ylim(0, 1)

(parsS, covs) = curve_fit(
    sigmoidDecay, 
    np.array([0]+distances), 
    np.array([1]+list(trps/sums))
)
(parsE, covs) = curve_fit(
    exponentialDecay, 
    np.array([0]+distances), 
    np.array([1]+list(trps/sums))
)
print(parsE)
samps = np.arange(0, 80, 1)
fitS = [sigmoidDecay(d, *parsS) for d in samps]
fitE = [exponentialDecay(d, *parsE) for d in samps]
(fig, ax) = plt.subplots(figsize=(10, 5))
plt.plot([0]+distances, [1]+list(trps/sums), 'o', label='data')
plt.plot(samps, fitS, '-', label='Sigmoid')
plt.plot(samps, fitE, '-', label='Exponential')
plt.legend()
ax.set_xlim(0, 80)
ax.set_ylim(0, 1)


###############################################################################
# Fit Data
###############################################################################
kernelDict = {
    'kernel': srv.exponentialDecay, 
    'params': {'A': 1, 'b': parsE[0]}
}
dHalf = srv.nSolveKernel(kernelDict, 0.5, guess=0, latlon=False)
