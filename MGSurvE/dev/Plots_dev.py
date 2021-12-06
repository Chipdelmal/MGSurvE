#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math
import numpy as np
import pandas as pd
from os import path
from sys import argv
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import MGSurvE as srv

(OUT_PTH, LND_TYPE, ID) = ('./Lands', 'UNIF', 'D01')
dta = srv.importLog(OUT_PTH, '{}_{}_LOG'.format(LND_TYPE, ID))
###############################################################################
# Plot
############################################################################### 
(fig, ax) = plt.subplots(figsize=(15, 15))
(fig, ax) = srv.plotGAEvolution(fig, ax, dta)
# srv.plotClean(fig, ax)
###############################################################################
# Save
############################################################################### 
pthSave = path.join(
    OUT_PTH, '{}_{}_GAP'.format(LND_TYPE, ID)
)
fig.savefig(
    pthSave, dpi=DPI, bbox_inches='tight', pad_inches=0, transparent=False
)
plt.close('all')


