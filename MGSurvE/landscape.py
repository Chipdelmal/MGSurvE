#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math
import numpy as np
import MGSurvE.matrices as mat
import MGSurvE.constants as cst
import MGSurvE.kernels as krn


class Landscape:
    """
    Stores the information for a mosquito landscape with different point-types 
        in the form of matrices and coordinates.

    Attributes
    ----------
    points : numpy array
        Sites coordinates in (x, y) or (lon, lat) format.
    pointTypes : numpy array
        Sites types in the same order and length as the points coordinates.
    pointsNumber: int
        Number of sites present in the environment.

    distanceMatrix : numpy array
        Matrix with the distances between all the points in the landscape.
    migrationMatrix : function
        Markov matrix that determines the probability of moving from one site
            to another.
    maskedMigrationMatrix : numpy array
        Markov matrix that biases migration probabilities as dictated by the 
            masking matrix.

    kernelFunction : function
        Function that determines de relationship between distances and 
            migration probabilities.
    kernelParams : dict
        Parameters required for the kernel function to determine migration
            probabilities.

    maskingMatrix : numpy array
        Matrix that determines the probability of shifting from one point-type
            to another one (squared with size equal to the number of point
            types)

    Methods
    -------
    calcPointsDistances()
        Calculates the distance matrix between the points in the landscape. Uses 
            the distanceFunction to calculate the distances. The attribute 
            distanceMatrix is updated in place.
    calcPointsMigration()
        Calculates the migration matrix between the points in the landscape 
            based on distance alone. The kernel function and kernel parameters
            attributes are used for the migration function. The attribute 
            migrationMatrix is updated in place.
    calcPointsMaskedMigration()
        Calculates the migration matrix that takes into account the point-types
            of the sites. Uses the masking matrix and migration matrix. The 
            attribute makedMigration is updated in place.
    """
    ###########################################################################
    # Initializers
    ###########################################################################
    def __init__(self, 
        points, 
        maskingMatrix=None,
        
        distanceMatrix=None, 
        distanceFunction=math.dist, 
        
        migrationMatrix=None, 
        kernelFunction=krn.zeroInflatedExponentialKernel,
        kernelParams={'params': cst.AEDES_EXP_PARAMS, 'zeroInflation': .75},

        maskedMigrationMatrix=None
    ):
        self.distanceFunction = distanceFunction
        self.kernelFunction = kernelFunction
        self.kernelParams = kernelParams
        self.maskingMatrix = maskingMatrix
        self.pointsNumber = len(points)
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
        if ('t' in ptsHead):
            self.pointTypes = np.asarray(points['t'])
        else:
            self.pointTypes = np.asarray([0]*len(points))
        # If no migration mask is provided, generate a dummy one --------------
        if maskingMatrix is None:
            ptNum = len(set(self.pointTypes))
            self.maskingMatrix = np.full((ptNum, ptNum), 1)
        else:
            self.maskingMatrix = np.asarray(maskingMatrix)
        # Init distance matrix ------------------------------------------------
        if distanceMatrix is None:
            self.calcPointsDistances()
        else:
            self.distanceMatrix = np.asarray(distanceMatrix)
        # Init migration matrix -----------------------------------------------
        if (migrationMatrix is None):
            self.calcPointsMigration()
        else:
            self.migrationMatrix = np.asarray(migrationMatrix)
        # Init masked migration matrix ----------------------------------------
        if (maskedMigrationMatrix is None) and (self.pointTypes is not None):
            self.calcPointsMaskedMigrationMatrix()
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

    def calcPointsMaskedMigrationMatrix(self):
         self.maskedMigration = mat.calcMaskedMigrationMatrix(
            self.migrationMatrix, self.maskingMatrix, self.pointTypes,
            distFun=self.distanceFunction
        )