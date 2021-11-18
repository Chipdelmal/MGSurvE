#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from deap import base, creator, algorithms, tools


def initChromosome(trapsNum, coordsRange):
    """ Generates a random uniform chromosome for GA optimization.
    
    Parameters:
        trapsNum (int): Number of traps to lay down in the landscape.
        xRan (tuple of tuples of floats): XY Range for the coordinates.

    Returns:
        (list): List of xy coordinates for the traps' positions.
    """
    (xRan, yRan) = coordsRange
    xCoords = np.random.uniform(xRan[0], xRan[1], trapsNum)
    yCoords = np.random.uniform(yRan[0], yRan[1], trapsNum)
    chromosome = [val for pair in zip(xCoords, yCoords) for val in pair]
    return chromosome


def reshapeInCanonicalForm(tau, sitesN, trapsN):
    """ Reshapes a migration matrix into canonical form (deprecated).
    
    Parameters:
        tau (numpy array): Traps migration matrix.
        sitesN (int): Number of sites.
        trapsN (int): Number of traps.

    Returns:
        (numpy array): Reshaped matrix in canonical form.
    """
    canO = list(range(sitesN, sitesN+trapsN))+list(range(0, sitesN))
    tauCan = np.asarray([[tau[i][j] for j in canO] for i in canO])
    return tauCan


def getMarkovAbsorbing(tauCan, trapsN):
    """ Get Markov's absorbing states (deprecated).
    
    Parameters:
        tauCan (numpy array): Traps migration matrix in canonical form.
        trapsN (int): Number of traps.

    Returns:
        (numpy array): Time to fall into absorbing states from anywhere in landscape.
    """
    A = tauCan[trapsN:, :trapsN]
    B = tauCan[trapsN:, trapsN:]
    F = np.linalg.inv(np.subtract(np.identity(B.shape[0]), B))
    return F


def getFundamentalMatrix(tau, sitesN, trapsN):
    """ Get Markov's fundamental matrix.
    
    Equivalent to using reshapeInCanonicalForm and getMarkovAbsorbing (which
        should be deprecated).

    Parameters:
        tau (numpy array): Traps migration matrix in canonical form.
        sitesN (int): Number of sites.
        trapsN (int): Number of traps.

    Returns:
        (numpy array): Time to fall into absorbing states from anywhere in landscape.
    """
    Q = tau[:sitesN, :sitesN]
    R = tau[:sitesN, -trapsN:]
    I = np.identity(Q.shape[0])
    F = np.linalg.inv(np.subtract(I, Q))
    return F