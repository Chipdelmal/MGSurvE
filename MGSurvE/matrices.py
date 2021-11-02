#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math
import numpy as np
import scipy.stats as stats
import MGSurvE.constants as cst
from sklearn.preprocessing import normalize


def calcDistanceMatrix(pointCoords, distFun=math.dist):
    coordsNum = len(pointCoords)
    distMatrix = np.empty((coordsNum, coordsNum))
    for (i, coordA) in enumerate(pointCoords):
        for (j, coordB) in enumerate(pointCoords):
            distMatrix[i][j] = distFun(coordA, coordB)
    return distMatrix


def calcMaskedMigration(
        migrationMatrix, maskingMatrix, pointTypes,
        distFun=math.dist
    ):
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