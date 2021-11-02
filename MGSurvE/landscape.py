#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math
import numpy as np
import MGSurvE.matrices as mat
import MGSurvE.constants as cst
import MGSurvE.kernels as krn


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
        kernelFunction=krn.zeroInflatedExponentialKernel,
        kernelParams={'params': cst.AEDES_EXP_PARAMS, 'zeroInflation': .75},

        maskedMigration=None
    ):
        self.distanceFunction = distanceFunction
        self.kernelFunction = kernelFunction
        self.kernelParams = kernelParams
        self.maskingMatrix = np.asarray(pointTypesMask)
        # Check and define coordinates ----------------------------------------
        ptsHead = set(points.columns)
        if ('x' in ptsHead) and ('y' in ptsHead):
            self.geometryType = 'xy'
            self.pointCoords = np.asarray(points[['x', 'y']])          
        elif ('lat' in ptsHead) and ('lon' in ptsHead):
            self.geometryType = 'll'
            self.pointCoords = np.asarray(points[['lon', 'lat']])
        else:
            raise Exception(
                ''' Check the landscape type! 
                Accepted headers are "(x, y)" and "(lat, lon).
                '''
            )
        # Check and define point-types ----------------------------------------
        self.pointTypes = None
        if ('t' in ptsHead):
            self.pointTypes = np.asarray(points['t'])
        # Init distance matrix ------------------------------------------------
        if distanceMatrix is None:
            self.calcPointsDistances()
        else:
            self.distanceMatrix = np.asarray(distanceMatrix)
        # Init migration matrix -----------------------------------------------
        if (migrationMatrix is None):
            self.calcPointsMigration()
        else:
            self.distanceMatrix = np.asarray(migrationMatrix)
        # Init masked migration matrix ----------------------------------------
        if (maskedMigration is None) and (self.pointTypes is not None):
            self.calcPointsMaskedMigration()
        else:
            self.maskedMigration = np.asarray(maskedMigration)
    ###########################################################################
    # Matrix Methods
    ###########################################################################
    def calcPointsDistances(self):
        self.distanceMatrix = mat.calcDistanceMatrix(
            self.pointCoords, self.distanceFunction
        )

    def calcPointsMigration(self):
        self.migrationMatrix = self.kernelFunction(
            self.distanceMatrix, **self.kernelParams
        )

    def calcPointsMaskedMigration(self):
         self.maskedMigration = mat.calcMaskedMigration(
            self.migrationMatrix, self.maskingMatrix, self.pointTypes,
            distFun=self.distanceFunction
        )