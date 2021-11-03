#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math
import numpy as np
import scipy.stats as stats
import MGSurvE.constants as cst
from sklearn.preprocessing import normalize

###############################################################################
# Migration Kernels
###############################################################################
def zeroInflatedLinearMigrationKernel(distMat, params=[.75, 1]):
    '''Calculates the zero-inflated linear distance-based movement probability.

    Args:
        distMat (numpy array): Distances matrix.

    Returns:
        numpy array: Migration matrix.
    '''
    coordsNum = len(distMat)
    migrMat = np.empty((coordsNum, coordsNum))
    for (i, row) in enumerate(distMat):
        for (j, dst) in enumerate(row):
            migrMat[i][j] = inverseLinearStep(dst, params=params)
        migrMat[i] = migrMat[i] / sum(migrMat[i])
    return migrMat


def truncatedExponential(distance, params=cst.AEDES_EXP_PARAMS):
    ''' Calculates the zero-inflated exponential distance-based movement probability.

    Args:
        distance (float): Distances matrix.
        params (list): Shape distribution parameters [rate, a, b].

    Returns:
        float: Migration probability.
    '''
    if(params[1] > params[2]):
        return None

    scale = 1.0/params[0]
    gA = stats.expon.cdf(params[1], scale=scale)
    gB = stats.expon.cdf(params[2], scale=scale)
    if np.isclose(gA, gB):
        return None

    densNum = stats.expon.pdf(distance, scale=scale)
    densDen = gB - gA

    return densNum/densDen


def zeroInflatedExponentialKernel(
        distMat, params=cst.AEDES_EXP_PARAMS, zeroInflation=.75
    ):
    '''Calculates the migration matrix using a zero-inflated exponential function.

    Args:
        distMat (numpy array): Distances matrix.
        params (list): Shape distribution parameters [rate, a, b].
        zeroInflation (float): Probability to stay in the same place.

    Returns:
        numpy array: Migration matrix.
    '''
    coordsNum = len(distMat)
    migrMat = np.empty((coordsNum, coordsNum))
    for (i, row) in enumerate(distMat):
        for (j, dst) in enumerate(row):
            if(i == j):
                migrMat[i][j] = 0
            else:
                migrMat[i][j] = truncatedExponential(dst, params=params)
        migrMat[i] = migrMat[i] / np.sum(migrMat[i]) * (1 - zeroInflation)
    np.fill_diagonal(migrMat, zeroInflation)
    tauN = normalize(migrMat, axis=1, norm='l1')
    return tauN