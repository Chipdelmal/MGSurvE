#!/usr/bin/env python
# -*- coding: utf-8 -*-

import subprocess
import numpy as np
import pandas as pd
from os import path
from sys import argv
from glob import glob
import matplotlib
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import MGSurvE as srv
import auxiliary as aux
from PIL import Image
matplotlib.use('agg')
plt.rcParams['axes.facecolor']='#00000000'
plt.rcParams['savefig.facecolor']='#00000000'


if srv.isNotebook():
    (ID, AP, RID) = ('YKNC', 'man', '02')
else:
    (ID, AP, RID) = argv[1:]
(GENS, RID, FPAT) = (5000, int(RID), ID+'-{}_16*')
###############################################################################
# File ID
###############################################################################
FNAME = FPAT.format(AP)[:-1]+'-{:02d}_'.format(RID)
I_PTH = '/home/chipdelmal/Documents/WorkSims/MGSurvE_Validations/YKND_5000'
O_PTH = path.join(I_PTH, 'VID_{}'.format(FNAME[:-1]))
srv.makeFolder(O_PTH)
###############################################################################
# Load Files
###############################################################################
(log, lnd) = (
    pd.read_csv(path.join(I_PTH, FNAME+'LOG.csv')),
    srv.loadLandscape(I_PTH, FNAME+'TRP', fExt='pkl')
)
###############################################################################
# Plot Clean Landscape
###############################################################################
(PROJ, FIGS, PAD, DPI) = (ccrs.PlateCarree(), (15, 15), 0, 300)
(fig, ax) = (plt.figure(figsize=FIGS), plt.axes(projection=PROJ))
lnd.plotSites(fig, ax, size=50)
srv.plotClean(fig, ax, bbox=lnd.landLimits)
fig.savefig(
    path.join(O_PTH, FNAME+'CLN.png'), 
    transparent=False, facecolor='w',
    bbox_inches='tight', pad_inches=PAD, dpi=DPI
)
plt.close("all")
###############################################################################
# Plot Optimization
###############################################################################
fitFun = aux.switchFunction(AP)
gen = 1000
for gen in range(GENS)[0:]:
    print("* Exporting {:04d}/{:04d}".format(gen, GENS), end='\r')
    # Get traps ---------------------------------------------------------------
    trpEntry = log.iloc[gen]['traps']
    if ID == 'YKND':
        trpPos = aux.idStringToArray(trpEntry, discrete=True)
        trpCds = srv.chromosomeIDtoXY(trpPos, lnd.pointID, lnd.pointCoords).T
    else:
        trpPos = aux.idStringToArray(trpEntry, discrete=False)
        trpCds = np.reshape(trpPos, (-1, 2)).T
    # Assemble and update traps -----------------------------------------------
    trapsLocs = pd.DataFrame(
        np.vstack([trpCds, lnd.trapsTypes, lnd.trapsFixed]).T, 
        columns=['lon', 'lat', 't', 'f']
    )
    trapsLocs['t'] = trapsLocs['t'].astype('int64')
    trapsLocs['f'] = trapsLocs['f'].astype('int64')
    lnd.updateTraps(trapsLocs, lnd.trapsKernels)
    lnd.updateTrapsRadii([0.250, 0.125, 0.100])
    # Calculate new fitness ---------------------------------------------------
    if ID == 'YKND':
        fitness = srv.calcDiscreteFitness(
            trpPos, lnd, optimFunction=srv.getDaysTillTrapped,
            optimFunctionArgs={'inner': np.sum, 'outer': fitFun}
        )[0]
    else:
        fitness = srv.calcFitness(
            trpPos, lnd, optimFunction=srv.getDaysTillTrapped,
            optimFunctionArgs={'inner': np.sum, 'outer': fitFun}
        )[0]
    # Plot --------------------------------------------------------------------
    (fig, ax) = (plt.figure(figsize=FIGS), plt.axes(projection=PROJ))
    # lnd.plotSites(fig, ax, size=50)
    lnd.plotTraps(
        fig, ax, 
        zorders=(30, 25), transparencyHex='55', 
        latlon=True, proj=PROJ
    )
    ax.text(
        0.8, 0.07, '{:05d}'.format(gen),
        horizontalalignment='right', verticalalignment='center',
        fontsize=25, color='#00000066',
        transform=ax.transAxes, zorder=5
    )
    ax.text(
        0.8, 0.10, '{:.02f}'.format(fitness),
        horizontalalignment='right', verticalalignment='center',
        fontsize=25, color='#00000066',
        transform=ax.transAxes, zorder=5
    )
    srv.plotClean(fig, ax, bbox=lnd.landLimits)
    fig.savefig(
        path.join(O_PTH, '{:04d}.png'.format(gen)), 
        transparent=True, facecolor=None,
        bbox_inches='tight', pad_inches=PAD, dpi=DPI
    )
    plt.close("all")
    # Overlay Brute-force -----------------------------------------------------
    # time.sleep(.1)
    background = Image.open(path.join(O_PTH, FNAME+'CLN.png')).convert('RGBA')
    foreground = Image.open(path.join(O_PTH, '{:04d}.png'.format(gen))).convert('RGBA')
    (w, h) = background.size
    background = background.crop((0, 0, w, h))
    foreground = foreground.resize((int(w/1), int(h/1)), Image.Resampling.LANCZOS)
    background = Image.alpha_composite(background, foreground)
    background.save(path.join(O_PTH, '{:04d}.png'.format(gen)), dpi=(DPI, DPI))
    background.close(); foreground.close()
###############################################################################
# Export FFMPEG
#   "ffmpeg -start_number 0 -r 4 -f image2 -s 1920x1080
#       -i STP_10_%05d.png -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" 
#       -vcodec libx264 -preset veryslow -crf 15 
#       -pix_fmt yuv420p OUTPUT_PATH.mp4"
############################################################################### 
fmpegBse = "ffmpeg -start_number 0 -r 12 -f image2 -s 1920x1080 -i {}/%04d.png ".format(O_PTH)
fmpegMid = "-vf pad=ceil(iw/2)*2:ceil(ih/2)*2 -pix_fmt yuv420p {}/{}.mp4 -y".format(O_PTH, FNAME)
fmpegFll = fmpegBse+fmpegMid
process = subprocess.Popen(fmpegFll.split(), stdout=subprocess.PIPE)
(output, error) = process.communicate()
    
