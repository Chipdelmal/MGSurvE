#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math
import numpy as np

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