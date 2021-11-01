#!/usr/bin/env python
# -*- coding: utf-8 -*-
###############################################################################
# Various auxiliary functions
###############################################################################

def makeFolder(path):
    """Crates a folder in the specified directory.

    Parameters
    ----------
    path : string
        Path of the folder than needs to be created.

    Returns
    -------
    NA
    """
    if not os.path.exists(path):
        try:
            os.mkdir(path)
        except OSError:
            raise OSError(
                    "Can't create destination directory (%s)!" % (path)
                )

def makeFolders(pathsList):
    """Crates a folders in the specified directories.

    Parameters
    ----------
    path : list
        List of paths of the folders than need to be created.

    Returns
    -------
    NA
    """
    for fldr in pathsList:
        makeFolder(fldr)


def isNotebook():
    """Checks if the script is running from a jupyter notebook.

    Parameters
    ----------
    NA

    Returns
    ----------
    NA
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