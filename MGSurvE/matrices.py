#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math
import numpy as np
import scipy.stats as stats
import MGSurvE.constants as cst
from sklearn.preprocessing import normalize


def calcDistanceMatrix(pointCoords, distFun=math.dist):
    """Calculates the distance matrix between all the provided coordinates.
    
    Args:
        pointCoords (numpy array): Coordinates of the sites.
        distFun (function): Distance function to be used in the computations.
    
    Returns:
        (numpy array): Distances matrix
    """
    coordsNum = len(pointCoords)
    distMatrix = np.empty((coordsNum, coordsNum))
    for (i, coordA) in enumerate(pointCoords):
        for (j, coordB) in enumerate(pointCoords):
            distMatrix[i][j] = distFun(coordA, coordB)
    return distMatrix


def calcMaskedMigrationMatrix(
        migrationMatrix, maskingMatrix, pointTypes
    ):
    """Calculates the masked migration matrix between points according to their types.
    
    Args:
        migrationMatrix (numpy array): Migration probabilities amongst points.
        maskingMatrix (numpy array): Transition probabilities between point-types.
        pointTypes (numpy vector): Point-types for each one of the sites in the matrix (in the same order).
    
    Returns:
        (numpy array): Masked migration matrix
    """
    pNum = len(migrationMatrix)
    itr = list(range(pNum))
    mskP = np.zeros((pNum, pNum))
    for row in itr:
        for col in itr:
            (a, b) = (pointTypes[row], pointTypes[col])
            mskP[row, col] = maskingMatrix[a, b]
    tauN = normalize(mskP*migrationMatrix, axis=1, norm='l1')
    return tauN


def calcTrapsToPointsDistances(trapsCoords, pointCoords, dFun=math.dist):
    """Generates the distances matrix between the traps and the sites.
    
    Args:
        trapsCoords (numpy array): Coordinates of the traps.
        pointsCoords (numpy array): Coordinates of the sites.
        dFun (function): Distance function to be used in the computations.
    
    Returns:
        (numpy array): Distances matrix
    """
    trapDists = np.asarray([
        [dFun(trap, site) for site in pointCoords]
        for trap in trapsCoords
    ]).T
    return trapDists


def calcTrapsProbabilities(
        trapsDistances, trapsTypes, trapsKernels, trapsMask, pointTypes
    ):
    """Calculates the traps probabilities given distances to points in the landscape and their effectiveness kernels.
    
    Args:
        trapsDistances (numpy array): Distances to all points.
        trapsTypes (numpy array): Types of the traps.
        trapsKernels (function): Kernels functions and params for trap types.
        trapsMask (np array): Traps' catching bias mask (with shape: trapTypes, pointTypes)
        pointTypes (list): 
    
    Returns:
        (numpy array): Traps probabilities
    """
    # Base unbiased probs -----------------------------------------------------
    trapProbs = np.asarray([
        [
            trapsKernels[ttype]['kernel'](i, **trapsKernels[ttype]['params']) 
            for (i, ttype) in zip(dist, trapsTypes)
        ] for dist in trapsDistances
    ])
    # Point-type to trap probs ------------------------------------------------
    trapMask = np.asarray([
        [
            trapsMask[ttype][pointTypes[ix]]
            for ttype in trapsTypes
        ] for ix in range(len(trapsDistances))
    ])
    # Alpha -------------------------------------------------------------------
    trapMatrix = trapProbs*trapMask
    return trapMatrix


def genVoidFullMigrationMatrix(migrationMatrix, trapsNumber):
    """Calculates a migration matrix with sections for traps (Xi) to be filled in place.
    
    Args:
        migrationMatrix (numpy array): Migration matrix without any traps (Tau).
        trapsNumber (int): Number of traps in the landscape
    
    Returns:
        (numpy array): Full migration matrix with no traps effects
    """
    (tau, sitesNumber) = (migrationMatrix, migrationMatrix.shape[0])
    identity = np.identity(trapsNumber)
    void = np.zeros((trapsNumber, sitesNumber))
    alpha = np.zeros((sitesNumber, trapsNumber))
    assembled = np.vstack([
        np.hstack([tau, alpha]), 
        np.hstack([void, identity])
    ])
    return assembled