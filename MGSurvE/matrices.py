#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math
import numpy as np
import scipy.stats as stats
import MGSurvE.constants as cst


###############################################################################
# Matrices and networks operations
###############################################################################
def calculateDistanceMatrix(landscape, distFun=math.dist):
    coordsNum = len(landscape)
    distMatrix = np.empty((coordsNum, coordsNum))
    for (i, coordA) in enumerate(landscape):
        for (j, coordB) in enumerate(landscape):
            distMatrix[i][j] = distFun(coordA, coordB)
    return distMatrix


###############################################################################
# Migration Kernels
###############################################################################
def zeroInflatedLinearMigrationKernel(
            distMat,
            params=[.75, 1]
        ):
    '''
    Takes in the distances matrix, zero inflated value (step) and two extra
        parameters to determine the change from distances into distance-based
        migration probabilities (based on the kernel function provided).
    '''
    coordsNum = len(distMat)
    migrMat = np.empty((coordsNum, coordsNum))
    for (i, row) in enumerate(distMat):
        for (j, dst) in enumerate(row):
            migrMat[i][j] = inverseLinearStep(dst, params=params)
        # Normalize rows to sum 1
        migrMat[i] = migrMat[i] / sum(migrMat[i])
    return migrMat


def truncatedExponential(distance, params=cst.AEDES_EXP_PARAMS):
    '''
    Calculates the zero-inflated exponential for the mosquito movement kernel
        (default parameters set to Aedes aegypti calibrations).
        params = [rate, a, b]
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
            distMat,
            params=cst.AEDES_EXP_PARAMS,
            zeroInflation=.75
        ):
    '''
    Calculates the migration matrix using a zero-inflated exponential function
        taking as arguments the species-specific lifespan parameters, and the
        kernel constants (along with the lifelong stay probability).
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
    return migrMat