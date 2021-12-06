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
x = range(dta.shape[0])    
(fig, ax) = plt.subplots(figsize=(15, 15))
plt.plot(x, dta['avg'], lw=2, color='#ffffffFF')
ax.fill_between(x, dta['min'], dta['max'], alpha=0.5, color='#1565c077', lw=0)
ax.set_xlim(0, max(x))
# ax.set_ylim(0, 5*minFits[-1])
ax.set_aspect((1/3)/ax.get_data_ratio())
###############################################################################
# Save
############################################################################### 
pthSave = path.join(
    PT_GA, '{}_{}-GA.png'.format(EXP_FNAME, str(TRAPS_NUM).zfill(2))
)
fig.savefig(
    pthSave, dpi=aux.DPI, bbox_inches='tight', pad_inches=0, transparent=False
)
plt.close('all')