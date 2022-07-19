'''Main object with sites, traps, and masks to optimize and visualize the landscape to be analyzed.

'''

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import dill
import math
import numpy as np
from time import time
from haversine import haversine
from sklearn.preprocessing import normalize
from sklearn.cluster import KMeans
import MGSurvE.matrices as mat
import MGSurvE.constants as cst
import MGSurvE.kernels as krn
import MGSurvE.plots as plt
import MGSurvE.optimization as opt
import MGSurvE.pointProcess as pts
import MGSurvE.auxiliary as aux
import warnings
warnings.filterwarnings('ignore', 'The iteration is not making good progress')

class Landscape:
    """ Stores the information for a mosquito landscape. Works with different point-types in the form of matrices and coordinates.
    
    Parameters:
        points (pandas dataframe): Sites positions with mandatory {x,y} coordinates and optional {t: type, a: attractiveness} values for each site in the landscape.

        kernelFunction (function): Function that determines de relationship between distances and migration probabilities.
        kernelParams (dict): Parameters required for the kernel function to determine migration probabilities.
        maskingMatrix (numpy array): Matrix that determines the probability of shifting from one point-type to another one (squared with size equal to the number of point types). If None, every point-type transition is equiprobable.
        attractionVector (list): Relative attraction of one site to another (does not preserve migration matrix diagonal)

        distanceMatrix (numpy array): Matrix with the distances between all the points in the landscape. If None, it's auto-calculated (see calcPointsDistances).
        migrationMatrix (numpy array): Markov matrix that determines the probability of moving from one site to another. If None, it's auto-calculated (see calcPointsMigration).
        maskedMigrationMatrix (numpy array): Markov matrix that biases migration probabilities as dictated by the masking matrix. If None, it's auto-calculated (see calcPointsMaskedMigration).

        distanceFunction (function): Function that takes two points in the landscape and calculates the distance between them.

        traps (pandas dataframe): Traps positions with mandatory {x,y} coordinates and optional {t: trap type integer, f: fixed bool (immovable)}.
        trapsKernels (dict): Traps' kernels functions and parameters in dictionary form (where the indices must match the "t" values in the traps dataframe).
        trapsMask (numpy array): Traps' masking matrix to make traps act upon mosquitos searching for a specific resource (shape: {trapsNum, sitesNum}).
        trapsRadii (list): List of probability values at which we want rings to be plotted in dataviz functions.

        landLimits (tuple): Landscape's bounding box.
    """
    ###########################################################################
    # Initializers
    ###########################################################################
    def __init__(self, 
        points,
        maskingMatrix=None,
        attractionVector=None,

        distanceMatrix=None, 
        distanceFunction=None, 
        
        migrationMatrix=None, 
        kernelFunction=krn.zeroInflatedExponentialKernel,
        kernelParams={'params': cst.SHORT_EXP_PARAMS, 'zeroInflation': .75},

        maskedMigrationMatrix=None,

        traps=None,
        trapsKernels={
            0: {
                'kernel': krn.exponentialDecay, 
                'params': cst.BASIC_EXP_TRAP,
                'inverse': None
            }
        },
        trapsMask=None,
        trapsRadii=[.25, .2, .1, .05],

        landLimits=None,
        populations=None
    ):
        """Constructor method
        """
        self.kernelFunction = kernelFunction
        self.kernelParams = kernelParams
        self.maskingMatrix = maskingMatrix
        self.attractionVector = attractionVector
        self.pointNumber = len(points)
        self.trapsCoords = None
        self.trapsTypes = None
        self.trapsDistances = None
        self.trapsKernels = trapsKernels
        self.trapsNumber = None
        self.trapsMigration = None
        self.trapsRadii = trapsRadii
        self.trapsFixed= None
        self.fundamentalMatrix = None
        self.trapsMask = None
        self.distanceFunction = distanceFunction
        self.populations = populations
        self.latlon = False
        self.landLimits = landLimits
        self.reversePlot = False
        # Check and define coordinates ----------------------------------------
        ptsHead = set(points.columns)
        if ('x' in ptsHead) and ('y' in ptsHead):
            self.geometryType = 'xy'
            self.pointCoords = np.asarray(points[['x', 'y']])
            if distanceFunction is None:
                self.distanceFunction = math.dist
        elif ('lat' in ptsHead) and ('lon' in ptsHead):
            self.geometryType = 'll'
            self.pointCoords = np.asarray(points[['lon', 'lat']])
            self.latlon = True
            self.reversePlot = True
            if distanceFunction is None:
                self.distanceFunction = lambda a, b: haversine(
                    (a[1], a[0]), (b[1], b[0]), unit='m'
                )
        else:
            raise Exception(
                '''Check the landscape type! 
                Accepted headers are (x, y) and (lon, lat).
                '''
            )
        # Check and define point-types ----------------------------------------
        if ('t' in ptsHead):
            self.pointTypes = np.asarray(points['t'])
        else:
            self.pointTypes = np.asarray([0]*len(points))
        # Check and define point-IDs ------------------------------------------
        if ('id' in ptsHead):
            self.pointID = list(points['id'])
        else:
            self.pointID = list(range(len(points)))
        # If no migration mask is provided, generate a dummy one --------------
        if maskingMatrix is None:
            ptNum = len(set(self.pointTypes))
            self.maskingMatrix = np.full((ptNum, ptNum), 1)
        else:
            self.maskingMatrix = np.asarray(maskingMatrix)
        # If no attraction vector is provided, generate a dummy one -----------
        if attractionVector is None:
            ptNum = self.pointNumber
            self.attractionVector = np.full(ptNum, 1)
        else:
            self.attractionVector = np.asarray(attractionVector)
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
            self.calcPointsMaskedMigration()
        else:
            self.maskedMigration = np.asarray(maskedMigrationMatrix)
        # Setup population sizes ----------------------------------------------
        if (self.populations is not None):
            if (self.populations.shape[0] != self.pointNumber):
                raise Exception(
                    '''Number of populations ({}) is not the same as 
                    number of sites ({}).
                    '''.format(self.populations.shape[0], self.pointNumber)
                )
        # Set land limits -----------------------------------------------------
        if self.landLimits is None:
            self.landLimits = self.getBoundingBox()
        else:
            self.landLimits = landLimits
        # Init traps locations ------------------------------------------------
        if (traps is not None):
            tpsHead = set(traps.columns)
            if (self.geometryType == 'xy'):
                self.trapsCoords = np.asarray(traps[['x', 'y']]).astype(float)
            else:
                self.trapsCoords = np.asarray(traps[['lon', 'lat']]).astype(float)
            # Check if there's trap-type information --------------------------
            if ('t' in tpsHead):
                self.trapsTypes = np.asarray(traps['t'])
            else:
                self.trapsTypes = np.asarray([0]*len(self.trapsCoords))
            if ('f' in tpsHead):
                self.trapsFixed = np.asarray(traps['f'])
            else:
                self.trapsFixed = np.asarray([0]*len(self.trapsCoords))
            if ('sid' in tpsHead):
                self.trapsSiteID = list(traps['sid'])
            if ('o' in tpsHead):
                self.trapsTOptim = np.asarray(traps['o'])
            else:
                self.trapsTOptim = np.asarray([0]*len(self.trapsCoords))
            self.trapsNumber = len(self.trapsCoords)
            # Init trap types -------------------------------------------------
            trapTypesNum = len(set(self.trapsTypes))
            pointTypesNum = len(set(self.pointTypes))
            if trapsMask is None:
                mskShape = (
                    max(set(self.trapsTypes))+1, 
                    len(set(self.pointTypes))
                )
                self.trapsMask = np.full(mskShape, 1)
            else:
                self.trapsMask = trapsMask
                if (self.trapsMask.shape[0] != trapTypesNum):
                    raise Exception(
                        '''Number of trap types in mask ({}) is not the same as 
                        trap types in traps dictionary ({}).
                        '''.format(self.trapsMask.shape[0], trapTypesNum)
                    )
                if (self.trapsMask.shape[1] != pointTypesNum):
                    raise Exception(
                        '''Number of point types in mask ({}) is not the same as 
                        the ones defined in the landscape ({}).
                        '''.format(self.trapsMask.shape[1], pointTypesNum)
                    )
            # Calculate trapsDistances ----------------------------------------
            self.calcTrapsDistances()
            # Generate empty traps matrix -------------------------------------
            self.trapsMigration = mat.genVoidFullMigrationMatrix(
                self.maskedMigration, self.trapsNumber
            )
            # Filll in traps matrix -------------------------------------------
            self.calcTrapsMigration()
            # Filll in traps matrix -------------------------------------------
            self.updateTrapsRadii(self.trapsRadii)
    ###########################################################################
    # Matrix Methods
    ###########################################################################
    def calcPointsDistances(self):
        """Calculates the distancesMatrix amongst the points (in place).
        """
        self.distanceMatrix = mat.calcDistanceMatrix(
            self.pointCoords, self.distanceFunction
        )

    def calcPointsMigration(self):
        """Calculates the migrationMatrix amongst the points (in place).
        """
        preAttractMigMatrix = self.kernelFunction(
            self.distanceMatrix, **self.kernelParams
        )
        self.migrationMatrix = mat.calcAttractiveness(
            preAttractMigMatrix, self.attractionVector
        )
    def calcPointsMaskedMigration(self):
        """Calculates the maskedMigrationMatrix depending on point-type (in place).
        """
        self.maskedMigration = mat.calcMaskedMigrationMatrix(
            self.migrationMatrix, self.maskingMatrix, self.pointTypes
        )
    ###########################################################################
    # Traps Methods
    ###########################################################################
    def calcTrapsDistances(self):
        """Calculates the trapsDistances matrix (in place).
        """
        self.trapsDistances = mat.calcTrapsToPointsDistances(
            self.trapsCoords, self.pointCoords, self.distanceFunction
        )
    def calcTrapsMigration(self):
        """Replaces section in the trapsMigration matrix (in place).
        """
        trapProbs = mat.calcTrapsProbabilities(
            self.trapsDistances, self.trapsTypes, 
            self.trapsKernels, self.trapsMask, self.pointTypes
        )
        ptsN = self.pointNumber
        self.trapsMigration[:ptsN, ptsN:] = trapProbs
        self.trapsMigration = normalize(self.trapsMigration, axis=1, norm='l1')
    def updateTrapsCoords(self, trapsCoords):
        self.trapsCoords = trapsCoords
        self.calcTrapsDistances()
        self.calcTrapsMigration()
    def updateTraps(self, traps, trapsKernels):
        """Updates the traps locations and migration matrices (in place).

        Parameters:
            traps (pandas dataframe):
            trapsKernel (dictionary):
        """
        oldTrapsNum = self.trapsNumber
        # Update traps locations and types ------------------------------------
        if (self.geometryType == 'xy'):
            self.trapsCoords = np.asarray(traps[['x', 'y']])
        else:
            self.trapsCoords = np.asarray(traps[['lon', 'lat']])
        self.trapsTypes = np.asarray(traps['t'])
        self.trapsNumber = len(self.trapsCoords)
        self.trapsKernels = trapsKernels
        self.trapsFixed = traps['f']
        try:
            self.trapsTOptim = traps['o']
        except:
            # print("No 'o' part of update")
            pass
        # Updating necessary matrices -----------------------------------------
        if (oldTrapsNum != len(self.trapsCoords)):
            self.trapsMigration = mat.genVoidFullMigrationMatrix(
                self.maskedMigration, self.trapsNumber
            )
            mskShape = (
                len(set(self.trapsTypes)), 
                len(set(self.pointTypes))
            )
            self.trapsMask = np.full(mskShape, 1)
        self.calcTrapsDistances()
        self.calcTrapsMigration()
        self.updateTrapsRadii(self.trapsRadii)
    def updateTrapsRadii(self, probValues):
        tker = self.trapsKernels
        for k in list(tker.keys()):
            tker[k].update({
                'radii': [
                    krn.nSolveKernel(tker[k], d, guess=0, latlon=self.latlon) 
                    for d in probValues
                ]
            })
        self.trapsKernels = tker
    ###########################################################################
    # Optimization
    ###########################################################################
    def calcFundamentalMatrix(self):
        """Calculates the Markov fundamental matrix on the landscape.
        """
        self.fundamentalMatrix = opt.getFundamentalMatrix(
            self.trapsMigration, self.pointNumber, self.trapsNumber
        )
    def getDaysTillTrapped(self, fitFuns={'outer': np.mean, 'inner': np.max}):
        """Gets the number of timesteps until a walker falls into a trap
        """
        if self.fundamentalMatrix is None:
            raise Exception(
                '''Calculate the Fundamental Matrix first! 
                lnd.calcFundamentalMatrix()
                '''
            )
        daysTillTrapped = opt.getFundamentalFitness(
            self.fundamentalMatrix, fitFuns
        )
        return daysTillTrapped
    ###########################################################################
    # Aggregate Landscape
    ###########################################################################
    # def aggregateLandscape(self, 
    #         clustersNumber, clusterAlgorithm=KMeans,
    #         randomState=time()
    #     ):
    #     clst = pts.clusterLandscape(
    #         self.pointCoords, clustersNumber,
    #         randomState=randomState, clusterAlgorithm=clusterAlgorithm
    #     )
    #     migMatAgg = pts.aggregateLandscape(
    #         self.migrationMatrix, clst['clusters']
    #     )
    #     self.pointCoords = clst['centroids']
    #     return migMatAgg
    ###########################################################################
    # Plotting Methods
    ###########################################################################
    def plotSites(self, 
            fig, ax, 
            markers=cst.MKRS, colors=cst.MCOL,
            size=350, edgecolors='w', linewidths=2,
            zorder=5, **kwargs
        ):
        """Plots the sites coordinates.
        """
        # pCoords = self.pointCoords
        # if reversePlot:
        #     pCoords = np.asarray([[i[1], i[0]] for i in self.pointCoords])
        (fig, ax) = plt.plotSites(
            fig, ax, 
            self.pointCoords, self.pointTypes,
            colors=colors, size=size, markers=markers,
            edgecolors=edgecolors, linewidths=linewidths,
            zorder=zorder, **kwargs
        )
        return (fig, ax)
    def plotMigrationNetwork(self,
            fig, ax, 
            lineColor='#03045e', lineWidth=20, 
            alphaMin=.5, alphaAmplitude=2.5,
            zorder=0, **kwargs
        ):
        """Plots the base sites migration network.
        """
        plt.plotMigrationNetwork(
            fig, ax, 
            self.migrationMatrix, self.pointCoords, self.pointCoords,
            lineColor=lineColor, lineWidth=lineWidth, 
            alphaMin=alphaMin, alphaAmplitude=alphaAmplitude,
            zorder=zorder, **kwargs
        )
        return (fig, ax)
    def plotMaskedMigrationNetwork(self,
            fig, ax, 
            lineColor='#03045e', lineWidth=20, 
            alphaMin=.5, alphaAmplitude=2.5,
            zorder=0, **kwargs
        ):
        """Plots the base sites migration network.
        """
        plt.plotMigrationNetwork(
            fig, ax, 
            self.maskedMigration, self.pointCoords, self.pointCoords,
            lineColor=lineColor, lineWidth=lineWidth, 
            alphaMin=alphaMin, alphaAmplitude=alphaAmplitude,
            zorder=zorder, **kwargs
        )
        return (fig, ax)
    def plotTraps(self,
            fig, ax, 
            colors=cst.TRP_COLS, marker="X",
            edgecolor='w', lws=(2, 0), ls=':',
            size=300, zorders=(25, -5),
            **kwargs
        ):
        """Plots the traps locations.
        """
        (fig, ax) = plt.plotTraps(
            fig, ax, 
            self.trapsCoords, self.trapsTypes, 
            self.trapsKernels, self.trapsFixed,
            colors=colors, marker=marker,
            edgecolor=edgecolor, lws=lws, ls=ls,
            size=size, zorders=zorders,
            **kwargs
        )
        return (fig, ax)
    def plotTrapsNetwork(self,
            fig, ax, 
            lineColor='#f72585', lineWidth=20, 
            alphaMin=.5, alphaAmplitude=1.5,
            zorder=0, **kwargs
        ):
        """Plots the traps networks.
        """
        (fig, ax) = plt.plotTrapsNetwork(
            fig, ax, 
            self.trapsMigration, self.trapsCoords, self.pointCoords,
            lineColor=lineColor, lineWidth=lineWidth, 
            alphaMin=alphaMin, alphaAmplitude=alphaAmplitude,
            zorder=zorder, **kwargs
        )
        return (fig, ax)
    def plotDirectedNetwork(self,
            fig, ax, 
            markers=cst.MKRS, colors=cst.MCOL,
            edgecolors='black'
        ):
        (fig, ax)= plt.plotDirectedNetwork(
            fig, ax, 
            sites=self.pointCoords, pTypes=self.pointTypes,
            transMtx=self.migrationMatrix,
            markers=markers, colors=colors, edgecolors=edgecolors
        )
        return (fig, ax)
    def plotLandBoundary(self,
            fig, ax, 
            landTuples=cst.LAND_TUPLES
        ):
        (fig, ax) = plt.plotLandBoundary(fig, ax, landTuples=landTuples)
        return (fig, ax)
    ###########################################################################
    # Auxiliary Methods
    ###########################################################################
    def getBoundingBox(self):
        """Returns the landscape's ((minX, maxX), (minY, maxY)).
        """
        (minX, minY) = np.apply_along_axis(min, 0, self.pointCoords)
        (maxX, maxY) = np.apply_along_axis(max, 0, self.pointCoords)
        retPair = ((minX, maxX), (minY, maxY))
        return retPair
