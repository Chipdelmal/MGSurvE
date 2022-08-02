'''GA Operators to calculate fitness and perform operations to search through optimization space.

'''

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import random
import numpy as np
import pandas as pd
from os import path
from random import choice
import numpy.random as rand
from deap import base, creator, algorithms, tools

###############################################################################
# Fitness function
###############################################################################
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
    """ Get Markov's fundamental matrix (pseudo-inverse).
    
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

def getFundamentalMatrixPseudoInverse(tau, sitesN, trapsN, rcond=1e-20):
    """ Get Markov's fundamental matrix (inverse).
    
    Equivalent to using reshapeInCanonicalForm and getMarkovAbsorbing (which
        should be deprecated).

    Parameters:
        tau (numpy array): Traps migration matrix in canonical form.
        sitesN (int): Number of sites.
        trapsN (int): Number of traps.
        rcond (float): Cutoff for small singular values.

    Returns:
        (numpy array): Time to fall into absorbing states from anywhere in landscape.
    """
    Q = tau[:sitesN, :sitesN]
    R = tau[:sitesN, -trapsN:]
    I = np.identity(Q.shape[0])
    F = np.linalg.pinv(np.subtract(I, Q), rcond=rcond)
    return F

def getFundamentalVector(tau, sitesN):
    # Equivalent to:
    #   np.sum(srv.getFundamentalMatrix(tau, sitesN, trapsN), axis=1)
    #   np.sum(srv.getFundamentalMatrixPseudoInverse(tau, sitesN, trapsN), axis=1)
    Q = tau[:sitesN, :sitesN]
    I = np.identity(Q.shape[0])
    o = np.ones(Q.shape[0])
    F = np.linalg.solve(I-Q, o)
    return F

def getFundamentalFitness(
        fundamentalMatrix, 
        fitFuns={'outer': np.mean, 'inner': np.max}
    ):
    """ Get fitness from Markov's fundamental matrix.

    Parameters:
        fundamentalMatrix (numpy array): Markov's fundamental matrix (calcFundamentalMatrix)
        fitFuns (dict): Dictionary containing the inner (row) and outer (col) operations for the fundamental matrix.

    Returns:
        (float): Summarized fitness function for the fundamental matrix.
    """
    if fundamentalMatrix.ndim == 2:
        # Using the Markov matrix
        daysInSites = np.apply_along_axis(
            fitFuns['inner'], 1, fundamentalMatrix
        )
    else:
        # Using the Markov vector
        daysInSites = fundamentalMatrix
    daysTillTrapped = fitFuns['outer'](daysInSites)
    return daysTillTrapped


###############################################################################
# GA (Basic)
###############################################################################
def initChromosome(trapsCoords, fixedTrapsMask, coordsRange):
    """ Generates a random uniform chromosome for GA optimization.
    
    Parameters:
        trapsCoords (int): Number of traps to lay down in the landscape.
        fixedTrapsMask (list of bools): Mask with coordinates that can be moved (true) and which can't (false).
        coordsRange (tuple of tuples of floats).
    Returns:
        (list): List of xy coordinates for the traps' positions.
    """
    (xRan, yRan) = coordsRange
    trapsNum = trapsCoords.shape[0]
    chromosome = trapsCoords.flatten()
    allele = 0
    for _ in range(trapsNum):
        if fixedTrapsMask[allele]:
            chromosome[allele+0] = np.random.uniform(xRan[0], xRan[1], 1)[0]
            chromosome[allele+1] = np.random.uniform(yRan[0], yRan[1], 1)[0]
        allele = allele + 2
    return chromosome

def initDiscreteChromosome(ptsIds, fixedTraps, trapsSiteID=None, banSites=None):
    """Generates a random uniform chromosome for discrete GA optimizations (from available sites).

    Args:
        ptsIds (list): List of points ids as stored in landscape (lnd.pointID).
        fixedTraps (list): List of flags on traps that can be moved (lnd.trapsFixed).
        banSites (set): List of sites that should not be taken into account for optimization (lnd.trapsBan).

    Returns:
        list: List of sites at which the traps are located for GA optimization.
    """    
    trapsNum = len(fixedTraps)  
    # If some of the sites are banned, update the list of valid points
    validNodes = ptsIds
    if banSites:
        validNodes = tuple(set(ptsIds)-banSites)
    # Iterate through the alleles to generate a random chromosome
    chromosome = [0]*trapsNum
    for ix in range(trapsNum):
        if not fixedTraps[ix]:
            chromosome[ix] = choice(validNodes)
        else:
            chromosome[ix] = trapsSiteID[ix] 
    return chromosome

def genFixedTrapsMask(trapsFixed, dims=2):
    """ Creates a mask for the fixed traps (non-movable).
    
    Parameters:
        trapsFixed (bool numpy array): Boolean array with the traps that are not movable (lnd.trapsFixed).
        dims (int): Unused for now, but it's the number of dimensions for the landscape.

    Returns:
        (numpy array): Mask of the elements that can be moved in the GA operations.
    """
    dups = [list([not(i)])*dims for i in trapsFixed]
    dupsVct = [item for sublist in dups for item in sublist]
    return np.asarray(dupsVct)

def mutateChromosome(
        chromosome, fixedTrapsMask,
        randFun=rand.normal, 
        randArgs={'loc': 0, 'scale': 0.1},
        indpb=0.5
    ):
    """ Mutates a chromosome with a probability distribution based on the mutation mask (in place).
    
    Parameters:
        chromosome (floats numpy array): GA's float chromosome generated by initChromosome.
        fxdTrpsMsk (bool numpy array): Array of bools that define which alleles can be mutated (1).
        randFun (function): Probability function for the mutation operation.
        randArgs (dict): Arguments to control the shape of the probability function.
        indpb (float): Independent probability to mutate each allele.

    Returns:
        (numpy array list): Selectively-mutated chromosome.
    """
    cLen = len(chromosome)
    randDraw = randFun(size=cLen, **randArgs)
    randMsk = randDraw*fixedTrapsMask
    for i in range(len(chromosome)):
        if (random.random() < indpb) and (fixedTrapsMask[i]):
            chromosome[i] = chromosome[i] + randMsk[i]
    return (chromosome, )

def mutateChromosomeAsymmetric(
        chromosome, fixedTrapsMask,
        randFun=rand.normal, 
        randArgs={
            'x': {'loc': 0, 'scale': 0.1}, 
            'y': {'loc': 0, 'scale': 0.1}
        },
        indpb=0.5
    ):
    """ Mutates a chromosome with a probability distribution based on the mutation mask with different probabilities for XY elements (in place).
    
    Parameters:
        chromosome (floats numpy array): GA's float chromosome generated by initChromosome.
        fxdTrpsMsk (bool numpy array): Array of bools that define which alleles can be mutated (1).
        randFun (function): Probability function for the mutation operation.
        randArgs (dict): Arguments to control the shape of the probability function ('x' and 'y' entries).
        indpb (float): Independent probability to mutate each allele.

    Returns:
        (numpy array list): Selectively-mutated chromosome.
    """
    cLen=len(chromosome)
    # Draw mutations for XY chromosomes ---------------------------------------
    randDrawX = randFun(size=int(cLen/2), **randArgs['x'])
    randDrawY = randFun(size=int(cLen/2), **randArgs['y'])
    # Interweave XY mutations -------------------------------------------------
    randDraw = np.empty((cLen,))
    randDraw[0::2] = randDrawX
    randDraw[1::2] = randDrawY
    # Mask and return results -------------------------------------------------
    randMsk = randDraw*fixedTrapsMask
    for i in range(len(chromosome)):
        if (random.random() < indpb) and (fixedTrapsMask[i]):
            chromosome[i] = chromosome[i] + randMsk[i]
    return (chromosome, )

def mutateDiscreteChromosome(
        chromosome, ptsIds, fixedTraps, 
        indpb=0.5, banSites=None
    ):
    """Mutates a discrete chromosome from the available sites in the landscape (in place).

    Args:
        chromosome (list): GA's float chromosome generated by initChromosome.
        ptsIds (list): List of IDs for the sites (points) in the landscape (lnd.pointID).
        fixedTraps (list): List of flags on traps that can be moved (lnd.trapsFixed).
        indpb (float, optional): Mutation probability for each independent allele. Defaults to 0.5.
        banSites (set, optional): List of sites that should not be taken into account for optimization (lnd.trapsBan). Defaults to None.

    Returns:
        list: Mutated chromosome (in place).
    """    
    # If some of the sites are banned, update the list of valid points
    validNodes = ptsIds
    if banSites:
        validNodes = tuple(set(ptsIds)-banSites)
    # Iterate through alleles for mutation
    cLen = len(fixedTraps)
    for i in range(cLen):
        if (random.random() < indpb) and not fixedTraps[i]:
            chromosome[i] = choice(validNodes)
    return (chromosome, )

def cxBlend(
        ind1, ind2, 
        fixedTrapsMask, 
        alpha=.5
    ):
    """ Mates two chromosomes by "blend" based on the provided mask (in place).
    
    This implementation is similar to DEAP's cxBlend (https://deap.readthedocs.io/en/master/api/tools.html#deap.tools.cxBlend). 
    Follow this link for the original code: https://github.com/DEAP/deap/blob/master/deap/tools/crossover.py

    Parameters:
        ind1 (floats numpy array): GA's float chromosome generated by initChromosome.
        ind2 (floats numpy array): GA's float chromosome generated by initChromosome.
        fxdTrpsMsk (bool numpy array): Array of bools that define which alleles can be mutated (1).
        alpha (float): weight for each of the chromosomes.

    Returns:
        (list of chromosomes): Mated individuals.
    """
    (offA, offB) = (ind1[:], ind2[:])
    for i, (x1, x2) in enumerate(zip(ind1, ind2)):
        if fixedTrapsMask[i]:
            gamma = (1. + 2. * alpha) * random.random() - alpha
            offA[i] = (1. - gamma) * x1 + gamma * x2
            offB[i] = gamma * x1 + (1. - gamma) * x2
    (ind1[:], ind2[:]) = (offA[:], offB[:])
    return (ind1, ind2)

def cxDiscreteUniform(ind1, ind2, fixedTraps, indpb=.5):
    """Mates two chromosomes in place by swapping alleles between them.

    Args:
        ind1 (list): Chromosome of the first parent in the mating operation.
        ind2 (list): Chromosome of the first parent in the mating operation.
        fixedTraps (_type_): List of flags on traps that can be moved (lnd.trapsFixed).
        indpb (float, optional): Mutation probability for each independent allele. Defaults to 0.5.

    Returns:
        (list of chromosomes): Mated individuals.
    """    
    (offA, offB) = (ind1[:], ind2[:])
    for ix in range(len(offA)):
        if not fixedTraps[ix] and (random.random() < indpb):
            offA[ix] = ind2[ix]
            offB[ix] = ind1[ix]
    (ind1[:], ind2[:]) = (offA[:], offB[:])
    return (ind1, ind2)
###############################################################################
# GA (Extended)
###############################################################################
def mutShuffleIndexes(individual, typeOptimMask, indpb=.5):
    (size, clen) = (len(typeOptimMask), len(individual))
    for i in range(size):
        # If the allele can be mutated and was sampled
        if (typeOptimMask[i]):
            if (random.random() < indpb):
                swap_indx = random.randint(0, clen-2)
                if swap_indx >= i:
                    swap_indx += 1
                # If sampIx is part of the extra pool, just swap
                if swap_indx >= size:
                    individual[i], individual[swap_indx] = individual[swap_indx], individual[i]
                # If sampIx is part of the placed traps, and it's optimizable
                elif typeOptimMask[swap_indx]:
                    individual[i], individual[swap_indx] = individual[swap_indx], individual[i]
                # If sampIx is part of the placed traps, but not optimizable
                else:
                    swap_indx = random.randint(size, clen-1)
                    individual[i], individual[swap_indx] = individual[swap_indx], individual[i]
    return (individual, )

def initChromosomeMixed(
        trapsCoords, 
        fixedTrapsMask, typeOptimMask,
        coordsRange, trapsPool, 
        indpb=.75
    ):
    coordSect = initChromosome(trapsCoords, fixedTrapsMask, coordsRange)
    typesInit = mutShuffleIndexes(trapsPool, typeOptimMask, indpb)[0]
    return [float(i) for i in coordSect]+list(typesInit)

def mutateChromosomeMixed(
        chromosome,
        fixedTrapsMask, typeOptimMask,
        mutCoordFun=mutateChromosome,
        mutCoordArgs={
            'randFun': rand.normal, 'randArgs': {'loc': 0, 'scale': 10}, 
            'indpb': 0.5
        },
        mutTypeFun=mutShuffleIndexes,
        mutTypeArgs={
            'indpb': 0.5
        }
    ):
    # Split chromosome in parts -----------------------------------------------
    (coordsSect, typesSect) = (
        chromosome[:len(fixedTrapsMask)], 
        chromosome[len(fixedTrapsMask):]
    )
    # Mutate coordinates section ----------------------------------------------
    coordsMut = mutCoordFun(coordsSect, fixedTrapsMask,**mutCoordArgs)[0]
    # Mutate types section ----------------------------------------------------
    typesMut = mutTypeFun(typesSect, typeOptimMask, **mutTypeArgs)[0]
    # Return mutated chromosome -----------------------------------------------
    return (coordsMut+typesMut)

###############################################################################
# Fitness Functions
###############################################################################
def getDaysTillTrapped(
        landscape,
        fitFuns={'outer': np.mean, 'inner': np.max}
    ):
    """Gets the number of timesteps until a walker falls into a trap.

    Parameters:
        landscape (object): Landscape object to use for the analysis.
        fitFuns (dict): Dictionary with the outer (row) and inner (col) functions to use on the matrix.

    Returns:
        (float): Number of days for mosquitoes to fall into traps given the fitFuns.
    """
    funMat = getFundamentalMatrix(
        landscape.trapsMigration, 
        landscape.pointNumber, 
        landscape.trapsNumber
    )   
    daysTillTrapped = getFundamentalFitness(funMat, fitFuns=fitFuns)
    return daysTillTrapped

def getDaysTillTrappedPseudoInverse(
        landscape,
        fitFuns={'outer': np.mean, 'inner': np.max},
        rcond=1e-30
    ):
    """Gets the number of timesteps until a walker falls into a trap (using pseudo-inverse matrix function).

    Parameters:
        landscape (object): Landscape object to use for the analysis.
        fitFuns (dict): Dictionary with the outer (row) and inner (col) functions to use on the matrix.
        rcond (float): Cutoff for small singular values.

    Returns:
        (float): Number of days for mosquitoes to fall into traps given the fitFuns.
    """
    funMat = getFundamentalMatrixPseudoInverse(
        landscape.trapsMigration, 
        landscape.pointNumber, 
        landscape.trapsNumber,
        rcond=rcond
    )     
    daysTillTrapped = getFundamentalFitness(funMat, fitFuns=fitFuns)
    return daysTillTrapped

def getDaysTillTrappedVector(
        landscape, fitFuns={'outer': np.mean, 'inner': None}
    ):
    funVct = getFundamentalVector(
        landscape.trapsMigration, 
        landscape.pointNumber
    )   
    daysTillTrapped = getFundamentalFitness(funVct, fitFuns=fitFuns)
    return daysTillTrapped

def calcFitness(
        chromosome, 
        landscape=None,
        optimFunction=getDaysTillTrapped,
        optimFunctionArgs={'outer': np.mean, 'inner': np.max}
    ):
    """Calculates the fitness function of the landscape given a chromosome (in place, so not thread-safe).

    Parameters:
        chromosome (floats numpy array): GA's float chromosome generated by initChromosome.
        landscape (object): Landscape object to use for the analysis.
        optimFunction (function): Function that turns a matrix into a fitness value.
        optimFunctionArgs (dict): Dictionary with the outer (row) and inner (col) functions to use on the matrix.

    Returns:
        (tuple of floats): Landscape's fitness function.
    """
    candidateTraps = np.reshape(chromosome, (-1, 2))
    landscape.updateTrapsCoords(candidateTraps)
    fit = optimFunction(landscape, fitFuns=optimFunctionArgs)
    return (float(abs(fit)), )

def calcFitnessPseudoInverse(
        chromosome, 
        landscape=None,
        optimFunction=getDaysTillTrappedPseudoInverse,
        optimFunctionArgs={'outer': np.mean, 'inner': np.max},
        rcond=1e-30
    ):
    """Calculates the fitness function of the landscape given a chromosome using the matrix pseudo-inverse function (in place, so not thread-safe).

    Parameters:
        chromosome (floats numpy array): GA's float chromosome generated by initChromosome.
        landscape (object): Landscape object to use for the analysis.
        optimFunction (function): Function that turns a matrix into a fitness value.
        optimFunctionArgs (dict): Dictionary with the outer (row) and inner (col) functions to use on the matrix.
        rcond (float): Cutoff for small singular values.

    Returns:
        (tuple of floats): Landscape's fitness function.
    """
    candidateTraps = np.reshape(chromosome, (-1, 2))
    landscape.updateTrapsCoords(candidateTraps)
    fit = optimFunction(landscape, fitFuns=optimFunctionArgs, rcond=rcond)
    return (float(abs(fit)), )

def calcSexFitness(
        chromosome, 
        landscapeMale=None, landscapeFemale=None,
        weightMale=1, weightFemale=1,
        optimFunction=getDaysTillTrapped,
        optimFunctionArgs={'outer': np.mean, 'inner': np.max}
    ):
    """Calculates the fitness function of a Male/Female set of landscapes with a weighted sum of the time-to catch between them.

    Parameters:
        chromosome (floats numpy array): GA's float chromosome generated by initChromosome.
        landscapeMale (object): Male landscape object to use for the analysis.
        landscapeFemale (object): Female landscape object to use for the analysis.
        weightMale (float): Preference on catching males over females.
        weightFemale (float): Preference on catching females over males.
        optimFunction (function): Function that turns a matrix into a fitness value.
        optimFunctionArgs (dict): Dictionary with the outer (row) and inner (col) functions to use on the matrix.

    Returns:
        (tuple of floats): Landscape's fitness function.
    """
    candidateTraps = np.reshape(chromosome, (-1, 2))
    landscapeMale.updateTrapsCoords(candidateTraps)
    landscapeFemale.updateTrapsCoords(candidateTraps)
    fit = [
        abs(optimFunction(lnd, fitFuns=optimFunctionArgs)) for lnd in 
        (landscapeMale, landscapeFemale)
    ]
    fitVal = (fit[0]*weightMale+fit[1]*weightFemale)/(2*(weightMale+weightFemale))
    return (fitVal, )

def chromosomeIDtoXY(chromosome, ptsID, pointCoords):
    siteIndex = [ptsID.index(i) for i in chromosome]
    trapXY = np.asarray([pointCoords[i] for i in siteIndex])
    return trapXY

def calcDiscreteFitness(
        chromosome, landscape,
        optimFunction=getDaysTillTrapped,
        optimFunctionArgs={'outer': np.mean, 'inner': np.max},
    ):
    """Calculates the fitness function of the landscape given a chromosome (in place, so not thread-safe).

    Parameters:
        chromosome (list): Discrete optimization chromosome.
        landscape (object): Landscape object to use for the analysis.
        optimFunction (function): Function that turns a matrix into a fitness value.
        optimFunctionArgs (dict): Dictionary with the outer (row) and inner (col) functions to use on the matrix.

    Returns:
        (tuple of floats): Landscape's fitness function.
    """
    trapXY = chromosomeIDtoXY(
        chromosome, landscape.pointID, landscape.pointCoords
    )
    fit = calcFitness(
        trapXY, landscape=landscape,
        optimFunction=optimFunction,
        optimFunctionArgs=optimFunctionArgs
    )
    return fit

def calcDiscreteFitnessPseudoInverse(
        chromosome, landscape,
        optimFunction=getDaysTillTrappedPseudoInverse,
        optimFunctionArgs={'outer': np.mean, 'inner': np.max},  
        rcond=1e-30
    ):
    """Calculates the fitness function of the landscape given a chromosome (in place, so not thread-safe).

    Parameters:
        chromosome (list): Discrete optimization chromosome.
        landscape (object): Landscape object to use for the analysis.
        optimFunction (function): Function that turns a matrix into a fitness value.
        optimFunctionArgs (dict): Dictionary with the outer (row) and inner (col) functions to use on the matrix.
        rcond (float): Cutoff for small singular values.

    Returns:
        (tuple of floats): Landscape's fitness function.
    """
    trapXY = chromosomeIDtoXY(
        chromosome, landscape.pointID, landscape.pointCoords
    )
    print(trapXY)
    fit = calcFitnessPseudoInverse(
        trapXY, landscape=landscape,
        optimFunction=optimFunction,
        optimFunctionArgs=optimFunctionArgs,
        rcond=rcond
    )
    return fit

def calcDiscreteSexFitness(
        chromosome, 
        landscapeMale=None, landscapeFemale=None,
        weightMale=1, weightFemale=1,
        optimFunction=getDaysTillTrapped,
        optimFunctionArgs={'outer': np.mean, 'inner': np.max}
    ):
    """Calculates the fitness function of a Male/Female set of landscapes with a weighted sum of the time-to catch between them.

    Parameters:
        chromosome (floats numpy array): GA's float chromosome generated by initChromosome.
        landscapeMale (object): Male landscape object to use for the analysis.
        landscapeFemale (object): Female landscape object to use for the analysis.
        weightMale (float): Preference on catching males over females.
        weightFemale (float): Preference on catching females over males.
        optimFunction (function): Function that turns a matrix into a fitness value.
        optimFunctionArgs (dict): Dictionary with the outer (row) and inner (col) functions to use on the matrix.

    Returns:
        (tuple of floats): Landscape's fitness function.
    """
    # Reshape traps positions into xy
    ptsIds = landscapeMale.pointID
    siteIndex = [ptsIds.index(i) for i in chromosome]
    trapXY = np.asarray([landscapeMale.pointCoords[i] for i in siteIndex])
    landscapeMale.updateTrapsCoords(trapXY)
    landscapeFemale.updateTrapsCoords(trapXY)
    fit = [
        abs(optimFunction(lnd, fitFuns=optimFunctionArgs)) for lnd in 
        (landscapeMale, landscapeFemale)
    ]
    fitVal = (fit[0]*weightMale+fit[1]*weightFemale)/(2*(weightMale+weightFemale))
    return (fitVal, )

###############################################################################
# Logging results
###############################################################################
def exportLog(
        logbook,
        outPath,
        filename
    ):
    """Dumps a dataframe with the report of the GA's history.

    Parameters:
        logbook (object): DEAP GA object.
        outPath (path): Path where the file will be exported.
        F_NAME (string): Filenamme (without extension).

    """
    if not isinstance(logbook, pd.DataFrame):
        log = pd.DataFrame(logbook)
    else:
        log = logbook
    log.to_csv(path.join(outPath, filename)+'.csv', index=False)

def importLog(
        inPath,
        filename
    ):
    """Gets the number of timesteps until a walker falls into a trap.

    Parameters:
        LOG_PTH (path): Path where the file is stored.
        F_NAME (dict): Filename with extension.

    Returns:
        (pandas dataframe): GA optimization log.
    """
    df = pd.read_csv(path.join(inPath, filename+'.csv'))
    return df

###############################################################################
# GA Wrapper
###############################################################################
def optimizeTrapsGA(
        landscape, 
        generations=1000,
        bbox='auto',pop_size='auto',
        mating_params={'mate': .3, 'cxpb': 0.5}, 
        mutation_params={'mean': 0, 'sd': 100, 'mutpb': .4, 'ipb': .5},
        selection_params={'tSize': 3},
        fitnessFun=calcFitness,
        optimFunction=getDaysTillTrapped, 
        fitFuns={'outer': np.mean, 'inner': np.max},
        verbose=True
    ):
    """Optimizes the traps' positions using a simple GA algorithm.

    Args:
        landscape (object): Landscape object to use for the analysis.
        generations (int, optional): Number of generations to run in the GA. Defaults to 1000.
        bbox (tuple, optional): If not 'auto', tuple with the landscape's bounding box for mutation operations. Defaults to 'auto'.
        pop_size (str, optional): If not 'auto', size of the chromosome population size in the GA. Defaults to 'auto'.
        mating_params (dict, optional): Mating probability ('mate') and crossover blending rate ('cxpb') for mating operations. Defaults to {'mate': .3, 'cxpb': 0.5}.
        mutation_params (dict, optional): Gaussian mean ('mean') and deviation ('sd') for mutation operations, as well as independent allele mutation probability ('ipb'). Defaults to {'mean': 0, 'sd': 100, 'mutpb': .4, 'ipb': .5}.
        selection_params (dict, optional): Tournament size for the selection algorithm. Defaults to {'tSize': 3}.
        optimFunction (function, optional): Fitness function to be used upon the movement matrices. Defaults to getDaysTillTrapped.
        fitFuns (dict, optional): Fitness matrix reduction statistics (inner applied first, and outter applied to the result). Defaults to {'outer': np.mean, 'inner': np.max}.
        verbose (bool, optional): Verbosity on the optimization. Defaults to True.

    Returns:
        (object, dataframe): Returns the landscape and logbook for the optimization.
    """    
    if bbox == 'auto':
        bbox = landscape.getBoundingBox()
    # GA parameters -----------------------------------------------------------
    if pop_size == 'auto':
        pop_size = int(10*(landscape.trapsNumber*1.25))
    if mating_params == 'auto':
        mating_params = {'mate': .3, 'cxpb': 0.5}
    if mutation_params == 'auto':
        mutation_params = {
            'mean': 0, 'sd': max([i[1]-i[0] for i in bbox])/2.5, 
            'mutpb': .4, 'ipb': .5
        }
    if selection_params == 'auto':
        selection_params = {'tSize': 3}
    trapsMask = genFixedTrapsMask(landscape.trapsFixed)
    ###########################################################################
    # Register GA Functions to DEAP
    ###########################################################################
    # Cost function to minimize -----------------------------------------------
    toolbox = base.Toolbox()
    creator.create("FitnessMin", base.Fitness, weights=(-1.0, ))
    creator.create("Individual", list, fitness=creator.FitnessMin)
    # Creators ----------------------------------------------------------------
    toolbox.register(
        "initChromosome", initChromosome, 
        trapsCoords=landscape.trapsCoords, 
        fixedTrapsMask=trapsMask, 
        coordsRange=bbox
    )
    toolbox.register(
        "individualCreator", tools.initIterate, 
        creator.Individual, toolbox.initChromosome
    )
    toolbox.register(
        "populationCreator", tools.initRepeat, 
        list, toolbox.individualCreator
    )
    # Mating and mutation operators -------------------------------------------
    toolbox.register(
        "mate", cxBlend, 
        fixedTrapsMask=trapsMask,
        alpha=mating_params['mate']
    )
    toolbox.register(
        "mutate", mutateChromosome,
        fixedTrapsMask=trapsMask,
        randArgs={'loc': mutation_params['mean'], 'scale': mutation_params['sd']}, 
        indpb=mutation_params['ipb']
    )
    # Select and evaluate -----------------------------------------------------
    toolbox.register(
        "select", tools.selTournament, 
        tournsize=selection_params['tSize']
    )
    toolbox.register(
        "evaluate", fitnessFun, 
        landscape=landscape,
        optimFunction=optimFunction,
        optimFunctionArgs=fitFuns
    )
    ###########################################################################
    # Registering GA stats
    ###########################################################################
    pop = toolbox.populationCreator(n=pop_size)
    hof = tools.HallOfFame(1)
    stats = tools.Statistics(lambda ind: ind.fitness.values)   
    stats.register("min", np.min)
    stats.register("avg", np.mean)
    stats.register("max", np.max)
    stats.register("std", np.std)
    stats.register(
        "best", lambda fitnessValues: fitnessValues.index(min(fitnessValues))
    )
    stats.register(
        "traps", lambda fitnessValues: pop[fitnessValues.index(min(fitnessValues))]
    )
    ###########################################################################
    # Optimization Cycle
    ###########################################################################
    (pop, logbook) = algorithms.eaSimple(
        pop, toolbox, ngen=generations,  
        cxpb=mating_params['cxpb'], mutpb=mutation_params['mutpb'],    
        stats=stats, halloffame=hof, verbose=verbose
    )
    ###############################################################################
    # Get and Export Results
    ############################################################################### 
    bestChromosome = hof[0]
    bestTraps = np.reshape(bestChromosome, (-1, 2))
    landscape.updateTrapsCoords(bestTraps)
    logDF = pd.DataFrame(logbook)
    return (landscape, logDF)


def optimizeTwoSexesTrapsGA(
        landscapeMale, landscapeFemale, sexWeights={'M': .5, 'F': .5},
        generations=1000,
        bbox='auto', pop_size='auto',
        mating_params={'mate': .3, 'cxpb': 0.5}, 
        mutation_params={'mean': 0, 'sd': 100, 'mutpb': .4, 'ipb': .5},
        selection_params={'tSize': 3},
        optimFunction=getDaysTillTrapped, 
        fitFuns={'outer': np.mean, 'inner': np.max},
        verbose=True
    ):
    """Optimizes the traps' positions using a simple GA algorithm for two-sexes kernels.

    Args:
        landscapeMale (object): Male landscape object to use for the analysis.
        landscapeFemale (object): Female landscape object to use for the analysis.
        sexWeights (dictionary): Male-to-Female priority dictionary.
        generations (int, optional): Number of generations to run in the GA. Defaults to 1000.
        bbox (tuple, optional): If not 'auto', tuple with the landscape's bounding box for mutation operations. Defaults to 'auto'.
        pop_size (str, optional): If not 'auto', size of the chromosome population size in the GA. Defaults to 'auto'.
        mating_params (dict, optional): Mating probability ('mate') and crossover blending rate ('cxpb') for mating operations. Defaults to {'mate': .3, 'cxpb': 0.5}.
        mutation_params (dict, optional): Gaussian mean ('mean') and deviation ('sd') for mutation operations, as well as independent allele mutation probability ('ipb'). Defaults to {'mean': 0, 'sd': 100, 'mutpb': .4, 'ipb': .5}.
        selection_params (dict, optional): Tournament size for the selection algorithm. Defaults to {'tSize': 3}.
        optimFunction (function, optional): Fitness function to be used upon the movement matrices. Defaults to getDaysTillTrapped.
        fitFuns (dict, optional): Fitness matrix reduction statistics (inner applied first, and outter applied to the result). Defaults to {'outer': np.mean, 'inner': np.max}.
        verbose (bool, optional): Verbosity on the optimization. Defaults to True.

    Returns:
        (object, dataframe): Returns the landscape and logbook for the optimization.
    """    
    if bbox == 'auto':
        bbox = landscapeMale.getBoundingBox()
    # GA parameters -----------------------------------------------------------
    if pop_size == 'auto':
        pop_size = int(10*(landscapeMale.trapsNumber*1.25))
    if mating_params == 'auto':
        mating_params = {'mate': .3, 'cxpb': 0.5}
    if mutation_params == 'auto':
        mutation_params = {
            'mean': 0, 'sd': max([i[1]-i[0] for i in bbox])/2.5, 
            'mutpb': .4, 'ipb': .5
        }
    if selection_params == 'auto':
        selection_params = {'tSize': 3}
    trapsMask = genFixedTrapsMask(landscapeMale.trapsFixed)
    ###########################################################################
    # Register GA Functions to DEAP
    ###########################################################################
    # Cost function to minimize -----------------------------------------------
    toolbox = base.Toolbox()
    creator.create("FitnessMin", base.Fitness, weights=(-1.0, ))
    creator.create("Individual", list, fitness=creator.FitnessMin)
    # Creators ----------------------------------------------------------------
    toolbox.register(
        "initChromosome", initChromosome, 
        trapsCoords=landscapeMale.trapsCoords, 
        fixedTrapsMask=trapsMask, 
        coordsRange=bbox
    )
    toolbox.register(
        "individualCreator", tools.initIterate, 
        creator.Individual, toolbox.initChromosome
    )
    toolbox.register(
        "populationCreator", tools.initRepeat, 
        list, toolbox.individualCreator
    )
    # Mating and mutation operators -------------------------------------------
    toolbox.register(
        "mate", cxBlend, 
        fixedTrapsMask=trapsMask,
        alpha=mating_params['mate']
    )
    toolbox.register(
        "mutate", mutateChromosome,
        fixedTrapsMask=trapsMask,
        randArgs={'loc': mutation_params['mean'], 'scale': mutation_params['sd']}, 
        indpb=mutation_params['ipb']
    )
    # Select and evaluate -----------------------------------------------------
    toolbox.register(
        "select", tools.selTournament, 
        tournsize=selection_params['tSize']
    )
    toolbox.register("evaluate", 
        calcSexFitness, 
        landscapeMale=landscapeMale,landscapeFemale=landscapeFemale,
        weightMale=sexWeights['M'], weightFemale=sexWeights['F'],
        optimFunction=optimFunction,
        optimFunctionArgs=fitFuns
    )
    ###########################################################################
    # Registering GA stats
    ###########################################################################
    pop = toolbox.populationCreator(n=pop_size)
    hof = tools.HallOfFame(1)
    stats = tools.Statistics(lambda ind: ind.fitness.values)   
    stats.register("min", np.min)
    stats.register("avg", np.mean)
    stats.register("max", np.max)
    stats.register("std", np.std)
    stats.register(
        "best", lambda fitnessValues: fitnessValues.index(min(fitnessValues))
    )
    stats.register(
        "traps", lambda fitnessValues: pop[fitnessValues.index(min(fitnessValues))]
    )
    ###########################################################################
    # Optimization Cycle
    ###########################################################################
    (pop, logbook) = algorithms.eaSimple(
        pop, toolbox, ngen=generations,  
        cxpb=mating_params['cxpb'], mutpb=mutation_params['mutpb'],    
        stats=stats, halloffame=hof, verbose=verbose
    )
    ###############################################################################
    # Get and Export Results
    ############################################################################### 
    bestChromosome = hof[0]
    bestTraps = np.reshape(bestChromosome, (-1, 2))
    landscapeMale.updateTrapsCoords(bestTraps)
    landscapeFemale.updateTrapsCoords(bestTraps)
    logDF = pd.DataFrame(logbook)
    return ((landscapeMale, landscapeFemale), logDF)
