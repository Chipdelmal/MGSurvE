#!/usr/bin/env python
# -*- coding: utf-8 -*-

import time
import numpy as np
import pandas as pd
import os
from os import path
from sys import argv
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from copy import deepcopy
import matplotlib.pyplot as plt
from compress_pickle import dump, load
import MGSurvE as srv
from PIL import Image

# ffmpeg -start_number START_NUMBER -r FRAMERATE -f image2 -s 1920x1080 -i STP_05_%05d.png -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -vcodec libx264 -preset veryslow -crf 15 -pix_fmt yuv420p OUTPUT_PATH.mp4


(OUT_PTH, LND_TYPE, ID) = (
    '/home/chipdelmal/Documents/WorkSims/MGSurvE_Benchmarks/STPVincenty/', # './Lands', 
    'STP', '01'
)
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
for i in range(0, gens):
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
    trapsLocs['t']=trapsLocs['t'].astype('int64')
    trapsLocs['f']=trapsLocs['f'].astype('int64')
    lnd.updateTraps(trapsLocs, lnd.trapsKernels)
    ###########################################################################
    # Plot Figure
    ###########################################################################
    (fig, ax) = plt.subplots(1, 1, figsize=(15, 15), sharey=False)
    lnd.plotTraps(fig, ax)
    lnd.plotLandBoundary(fig, ax)
    srv.plotClean(fig, ax, bbox=lnd.landLimits)
    ax.text(
        0.5, 0.5, '{:.3f}'.format(gaMin[i]),
        horizontalalignment='center', verticalalignment='center',
        fontsize=100, color='#00000011',
        transform=ax.transAxes, zorder=5
    )
    ax.text(
        0.5, 0.475, '{:03d}'.format(i),
        horizontalalignment='center', verticalalignment='center',
        fontsize=25, color='#00000011',
        transform=ax.transAxes, zorder=5
    )
    pthSave = path.join(IMG_PTH, '{}{:05d}.png'.format(fPat, i))
    fig.savefig(
        pthSave, 
        dpi=DPI, bbox_inches='tight', 
        pad_inches=0.1, transparent=True
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
    plt.close('all')
