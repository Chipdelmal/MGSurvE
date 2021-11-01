#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math
import numpy as np
import MGSurvE.matrices as mat


class Landscape:
    ###########################################################################
    # Default parameters
    ###########################################################################

    ###########################################################################
    # Initializers
    ###########################################################################
    def __init__(self, 
        pointCoordinates, pointTypes=None, pointTypesMask=None,
        trapsCoordinates=None, trapsTypes=None,
        distanceMatrix=None, migrationMatrix=None, 
        distanceFunction=math.dist, geometryType='xy'
    ):
        (self.distanceFunction, self.geometryType) = (
            distanceFunction, geometryType
        )
        # Points --------------------------------------------------------------
        # Check if the number of points is equal to the number of point-types
        #   provided or generate a default type (0) for them.
        pcNum = len(pointCoordinates)
        if pointTypes is None:
            pointTypes = [0] * pcNum
            pointTypesMask = np.asarray([[1]])
        (ptNum, ptmNum) = (len(pointTypes), len(set(pointTypes)))
        # If there's a missmatch in the lengths, raise an exception
        if (pcNum == ptNum) and (len(pointTypesMask) == ptmNum):
            self.pointCoordinates = pointCoordinates
            self.pointTypes = pointTypes
            self.pointTypesMask = pointTypesMask
        else:
            msg = 'Number of points ({}) does not match number of point types ({}) or point-type mask ({}).'
            raise Exception(msg.format(pcNum, ptNum, ptmNum))
        # Traps ---------------------------------------------------------------
        if trapsCoordinates is not None:
            (tcNum, ttNum) = (len(trapsCoordinates), len(trapsTypes))
            if pcNum == ptNum:
                self.pointCoordinates = pointCoordinates
                self.pointTypes = pointTypes
            else:
                msg = 'Number of traps ({}) does not match number of trap types ({}).'
                raise Exception(msg.format(tcNum, ttNum))
        # Matrices ------------------------------------------------------------
        if distanceMatrix is None:
            self.calculatePointsDistances()

    ###########################################################################
    # Matrix Methods
    ###########################################################################
    def calculatePointsDistances(self):
        self.distanceMatrix = mat.calculateDistanceMatrix(
            self.pointCoordinates, self.distanceFunction
        )