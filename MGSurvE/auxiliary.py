#!/usr/bin/env python
# -*- coding: utf-8 -*-

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