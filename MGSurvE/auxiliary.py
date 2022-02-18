
'''Various functions that are serve I/O purposes or not relevant to any specific module.

'''

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import numpy as np
import pandas as pd
from os import path
from compress_pickle import dump, load
from vincenty import vincenty


###############################################################################
# Vincenty distance between points
###############################################################################
def vincentyDistance(pointA, pointB, meters=False):
    """Calculates the Vincenty arc distance between points.
    
    Args:
        pointA (tuple): First point in (lon, lat) format.
        pointB (tuple): Second point in (lon, lat) format.

    Returns:
        (float): Distance between points
    """
    distKM = vincenty(
        (pointA[1], pointA[0]), 
        (pointB[1], pointB[0])
    )
    if meters:
        return distKM*1000
    else:
        return distKM


###############################################################################
# Various auxiliary functions
###############################################################################
def makeFolder(path):
    """Creates a folder if it doesn't already exist.
    
    Args:
        path (string): Path to the desired directory.
    """
    if not os.path.exists(path):
        try:
            os.mkdir(path)
        except OSError:
            raise OSError(
                    "Can't create destination directory (%s)!" % (path)
                )

def makeFolders(pathsList):
    """Creates a list of folders if they don't exist.
    
    Args:
        paths (list): List path to the desired directories.
    """
    for fldr in pathsList:
        makeFolder(fldr)


def isNotebook():
    """Checks if the script is running from a Jupyter environment.

    Returns:
        bool: Flags Jupyter environment. 
    """
    try:
        shell = get_ipython().__class__.__name__
        if shell == 'ZMQInteractiveShell':
            return True   # Jupyter notebook or qtconsole
        elif shell == 'TerminalInteractiveShell':
            return False  # Terminal running IPython
        else:
            return False  # Other type (?)
    except NameError:
        return False      # Probably standard Python interpreter

###############################################################################
# Dump/Load functions
###############################################################################
def dumpLandscape(landscape, fPath, fName, fExt='bz2'):
    """Exports a serialized landscape to disk.

    Args:
        landscape (object): Landscape object to export.
        fPath (path): Path to where the landscape will be exported.
        fName (string): Filename.
        fExt (string): File extension.
    """
    dump(
        landscape, 
        path.join(fPath, '{}.{}'.format(fName, fExt))
    )

def loadLandscape(fPath, fName, fExt='bz2'):
    """Loads a serialized landscape from disk.

    Args:
        fPath (path): Path from where the landscape will be loaded.
        fName (string): Filename.
        fExt (string): File extension.

    Returns:
        (object): Landscape object.
    """
    lnd = load(
        path.join(fPath, '{}.{}'.format(fName, fExt))
    )
    return lnd

def exportLandscape(landscape, fPath, fName):
    """Exports a landscape to disk to CSV files (for use in MGDrivE).

    Args:
        landscape (object): Landscape object to export.
        fPath (path): Path to where the landscape will be exported.
        fName (string): Filename.

    """
    coords = np.concatenate(
        (landscape.pointCoords, landscape.trapsCoords), 
        axis=0
    )
    sitesTypes = np.concatenate(
        (landscape.pointTypes, np.asarray(landscape.trapsNumber*[-1])),
        axis=0
    )
    trapTypes = np.concatenate(
        (np.asarray(landscape.pointNumber*[-1]), landscape.trapsTypes),
        axis=0
    )
    coordsT = coords.T
    outArray = np.vstack((coordsT[0], coordsT[1], sitesTypes, trapTypes)).T
    outDF = pd.DataFrame(
        outArray, columns=('x', 'y', 'SitesType', 'TrapsType')
    )
    # Export to disk
    outDF.to_csv(
        path.join(fPath, '{}_Sites.csv'.format(fName)),
        index=False
    )
    np.savetxt(
        path.join(fPath, '{}_Migration.csv'.format(fName)),
        landscape.trapsMigration, 
        delimiter=","
    )

###############################################################################
# Paths functions
###############################################################################
def makeFolder(path):
    """Crates a folder in the specified directory.

    Args:
        path (string): Path of the folder than needs to be created.

    """
    if not os.path.exists(path):
        try:
            os.mkdir(path)
        except OSError:
            raise OSError(
                    "Can't create destination directory (%s)!" % (path)
                )