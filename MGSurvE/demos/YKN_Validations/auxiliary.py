#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

def switchFunction(fid):
    if fid=='sum':
        return np.sum
    elif fid=='man':
        return np.mean
    else:
        return np.max
   
    
def idStringToArray(string):
    return np.array([int(i) for i in string[1:-1].split(',')])


def getBestTraps(log):
    idsStr = log.sort_values('min').iloc[0]
    (fit, trps) = (idsStr['min'], idsStr['traps'])
    trpsArray = idStringToArray(trps)
    return (fit, trpsArray)