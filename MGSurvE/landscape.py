#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math
import numpy as np
import MGSurvE.matrices as mat
import MGSurvE.constants as cst
import MGSurvE.kernels as krn


class Landscape:
    """ Stores the information for a mosquito landscape. Works with different point-types in the form of matrices and coordinates.
    
    Parameters:
        points (pandas dataframe): Sites coordinates in (x, y) or (lon, lat) format.
        pointTypes (numpy array): Sites types in the same order and length as the points coordinates.
               
        kernelFunction (function): Function that determines de relationship between distances and migration probabilities.
        kernelParams (dict): Parameters required for the kernel function to determine migration probabilities.
        maskingMatrix (numpy array): Matrix that determines the probability of shifting from one point-type to another one (squared with size equal to the number of point types). If None, every point-type transition is equiprobable.

        distanceMatrix (numpy array): Matrix with the distances between all the points in the landscape. If None, it's auto-calculated (see calcPointsDistances).
        migrationMatrix (numpy array): Markov matrix that determines the probability of moving from one site to another. If None, it's auto-calculated (see calcPointsMigration).
        maskedMigrationMatrix (numpy array): Markov matrix that biases migration probabilities as dictated by the masking matrix. If None, it's auto-calculated (see calcPointsMaskedMigration).

    Attributes:
        pointsNumber (int): Number of sites present in the environment.

        distanceMatrix (numpy array): Distances amongst the points in the landscape.
        migrationMatrix (numpy array): Distance-based migration probabilities amongst the points in the landscape.
        maskedMigrationMatrix (numpy array): Point-type based migration probabilities amongst the points in the landscape.

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
        """Constructor method
        """
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
        """Calculates the distancesMatrix amongst the points (in place). Uses the distanceFunction to calculate the distanceMatrix internally.
        """
        self.distanceMatrix = mat.calcDistanceMatrix(
            self.pointCoords, self.distanceFunction
        )

    def calcPointsMigration(self):
        """Calculates the migrationMatrix amongst the points (in place). Uses the kernelFunction and kernelParams to generate the migrationMatrix.
        """
        self.migrationMatrix = self.kernelFunction(
            self.distanceMatrix, **self.kernelParams
        )

    def calcPointsMaskedMigration(self):
        """Calculates the maskedMigrationMatrix depending on point-type (in place). Uses the maskingMatrix to bias the migrationMatrix and store the results in the maskedMigrationMatrix.
        """
        self.maskedMigration = mat.calcMaskedMigrationMatrix(
            self.migrationMatrix, self.maskingMatrix, self.pointTypes,
            distFun=self.distanceFunction
        )