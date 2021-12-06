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

(OUT_PTH, LND_TYPE, ID) = ('./Lands', 'UNIF', 'D01')
dta = srv.importLog(OUT_PTH, '{}_{}_LOG'.format(LND_TYPE, ID))