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
    print(migrationMatrix)
    for row in itr:
        for col in itr:
            (a, b) = (pointTypes[row], pointTypes[col])
            mskP[row, col] = maskingMatrix[a, b]
    tauN = normalize(mskP*migrationMatrix, axis=1, norm='l1')
    return tauN