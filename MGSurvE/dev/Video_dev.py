#!/usr/bin/env python
# -*- coding: utf-8 -*-

import time
import numpy as np
import pandas as pd
import os
from os import path
from sys import argv
from copy import deepcopy
import matplotlib.pyplot as plt
from compress_pickle import dump, load
import MGSurvE as srv
from PIL import Image


(OUT_PTH, LND_TYPE, ID) = ('./Lands', 'UNIF', '0D1')
fPat = '{}_{}_'.format(LND_TYPE, ID)
IMG_PTH = path.join(OUT_PTH, fPat+'VID')
srv.makeFolder(IMG_PTH)
DPI = 300

lnd = srv.loadLandscape(OUT_PTH, fPat+'TRP')
dat = srv.importLog(OUT_PTH, fPat+'LOG')
###############################################################################
# Plot Loop
############################################################################### 
(gaMin, gaTraps, gens) = (dat['min'], dat['traps'], dat.shape[0])
bbox = lnd.getBoundingBox()
i=10
for i in range(gens):
    print("* Exporting frame {:05d}".format(i), end='\r')
    ###########################################################################
    # Reshape and update traps
    ###########################################################################
    trapsCoords = np.reshape(
        np.fromstring(gaTraps[i][1:-1], sep=','), (-1, 2)
    ).T
    trapsLocs = pd.DataFrame(
        np.vstack([trapsCoords, lnd.trapsTypes, lnd.trapsFixed]).T, 
        columns=['x', 'y', 't', 'f']
    )
    lndUp = lnd.updateTraps(trapsLocs, lnd.trapsKernels)
    ###########################################################################
    # Plot Figure
    ###########################################################################
    (fig, ax) = plt.subplots(1, 1, figsize=(15, 15), sharey=False)
    lnd.plotTraps(fig, ax)
    srv.plotClean(fig, ax, frame=False)
    ax.set_xlim(*bbox[0])
    ax.set_ylim(*bbox[1])
    ax.text(
        0.5, 0.5, '{:.3f}'.format(gaMin[i]),
        horizontalalignment='center', verticalalignment='center',
        fontsize=100, color='#00000011',
        transform=ax.transAxes, zorder=5
    )
    pthSave = path.join(IMG_PTH, '{}{:05d}.png'.format(fPat, i))
    fig.savefig(
        pthSave, 
        dpi=DPI, bbox_inches='tight', 
        pad_inches=.1, transparent=True
    )
    plt.close('all')
    ###########################################################################
    # Overlay Brute-force
    ###########################################################################
    time.sleep(2)
    background = Image.open(path.join(OUT_PTH, fPat+'CLN.png'))
    foreground = Image.open(pthSave)
    (w, h) = background.size
    background = background.crop((0, 0, w, h))
    foreground = foreground.resize((int(w/1), int(h/1)),Image.ANTIALIAS)
    background.paste(foreground, (0, 0), foreground)
    background.save(pthSave, dpi=(DPI, DPI))
    background.close()
    foreground.close()
