#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math
import numpy as np
import MGSurvE.matrices as mat
import MGSurvE.constants as cst
 

class Landscape:
    ###########################################################################
    # Initializers
    ###########################################################################
    def __init__(self, 
        points, 
        pointTypesMask=None,
        
        distanceMatrix=None, 
        distanceFunction=math.dist, 
        
        migrationMatrix=None, 
        kernelFunction=mat.zeroInflatedExponentialKernel,
        kernelParams={'params': cst.AEDES_EXP_PARAMS, 'zeroInflation': .75}
    ):
        self.distanceFunction = distanceFunction
        self.kernelFunction = kernelFunction
        self.kernelParams = kernelParams
        # Check landscape type and define coordinates -------------------------
        ptsHead = set(points.columns)
        if ('x' in ptsHead) and ('y' in ptsHead) and ('t' in ptsHead):
            self.geometryType = 'xy'
            self.pointCoords = np.asarray(points[['x', 'y']])
            self.pointTypes = np.asarray(points['t'])
        elif ('lat' in ptsHead) and ('lon' in ptsHead) and ('t' in ptsHead):
            self.geometryType = 'll'
            self.pointCoords = np.asarray(points[['lon', 'lat']])
            self.pointTypes = np.asarray(points['t'])
        else:
            raise Exception(
                ''' Check the landscape type! 
                Accepted headers are "(x, y, t)" and "(lat, lon, t).
                '''
            )
        # Init distance matrix ------------------------------------------------
        if distanceMatrix is None:
            self.calculatePointsDistances()
        else:
            self.distanceMatrix = distanceMatrix
        # Init migration matrix -----------------------------------------------
        if migrationMatrix is None:
            self.calculatePointsMigration()
        else:
            self.distanceMatrix = distanceMatrix
    ###########################################################################
    # Matrix Methods
    ###########################################################################
    def calculatePointsDistances(self):
        self.distanceMatrix = mat.calculateDistanceMatrix(
            self.pointCoords, self.distanceFunction
        )

    def calculatePointsMigration(self):
        self.migrationMatrix = self.kernelFunction(
            self.distanceMatrix, **self.kernelParams
        )