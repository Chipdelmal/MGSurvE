#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math
import numpy as np
import pandas as pd
from os import path
from sys import argv
from copy import deepcopy
import matplotlib.pyplot as plt
from deap import base, creator, algorithms, tools
from compress_pickle import dump, load
import MGSurvE as srv
import warnings
warnings.filterwarnings('ignore', 'The iteration is not making good progress')


(OUT_PTH, LND_TYPE, ID) = ('./Lands', 'AGG', '001')

ptsNum = 400
bbox = ((-225, 225), (-175, 175))
xy = srv.ptsRandUniform(ptsNum, bbox).T
points = pd.DataFrame({'x': xy[0], 'y': xy[1], 't': [0]*xy.shape[1]})
lnd = srv.Landscape(
    points, 
    kernelParams={'params': [.075, 1.0e-10, math.inf], 'zeroInflation': .75}
)
aggLnd = deepcopy(lnd)
lnd.aggregateLandscape(100).shape