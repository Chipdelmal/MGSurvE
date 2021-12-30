#!/usr/bin/env python
# -*- coding: utf-8 -*-

from sys import argv
import math
from vincenty import vincenty
import numpy as np
import pandas as pd
from os import path
from sys import argv
from copy import deepcopy
import matplotlib.pyplot as plt
from deap import base, creator, algorithms, tools
from compress_pickle import dump, load
from sklearn.preprocessing import normalize
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import MGSurvE as srv
import warnings
warnings.filterwarnings('ignore', 'The iteration is not making good progress')


(ID, OUT_PTH) = (
    'STP', 
    # '/RAID5/marshallShare/MGS_Benchmarks/STP/'
    '/home/chipdelmal/Documents/WorkSims/MGSurvE_Benchmarks/STP/'
)
TRPS_NUM = 5# int(argv[1]) # 3
IX_SPLIT = 27
DIAG_VAL = 0 # 0.05
###############################################################################
# Load Pointset
###############################################################################
sites = pd.read_csv(path.join(OUT_PTH, 'stp_cluster_sites_pop_v5_fixed.csv'))
sites['t'] = [0]*sites.shape[0]
SAO_TOME_LL = sites.iloc[IX_SPLIT:]
SAO_bbox = (
    (min(SAO_TOME_LL['lon']), max(SAO_TOME_LL['lon'])),
    (min(SAO_TOME_LL['lat']), max(SAO_TOME_LL['lat']))
)
SAO_TOME_LL = SAO_TOME_LL .rename(
    columns={'lon': 'x', 'lat': 'y'}
)
###############################################################################
# Load Migration Matrix
###############################################################################
migration = np.genfromtxt(
    path.join(OUT_PTH, 'kernel_cluster_v6a.csv'), delimiter=','
)
msplit = migration[IX_SPLIT:,IX_SPLIT:]
np.fill_diagonal(msplit, DIAG_VAL)
SAO_TOME_MIG = normalize(msplit, axis=1, norm='l1')
###############################################################################
# Defining Traps
###############################################################################
nullTraps = [0]*TRPS_NUM
(lonTrap, latTrap) = (
    np.random.uniform(SAO_bbox[0][0], SAO_bbox[0][1], TRPS_NUM),
    np.random.uniform(SAO_bbox[1][0], SAO_bbox[1][1], TRPS_NUM)
)
traps = pd.DataFrame({
    'x': lonTrap, 'y': latTrap,
    't': nullTraps, 'f': nullTraps
})
tKer = {0: {'kernel': srv.exponentialDecay, 'params': {'A': .5, 'b': 100}}}
###############################################################################
# Setting Landscape Up
###############################################################################
lnd = srv.Landscape(
    SAO_TOME_LL, migrationMatrix=SAO_TOME_MIG,
    traps=traps, trapsKernels=tKer,
    distanceFunction=math.dist
)
bbox = lnd.getBoundingBox()
trpMsk = srv.genFixedTrapsMask(lnd.trapsFixed)

###############################################################################
# Plot
############################################################################### 
(fig, ax) = (
    plt.figure(figsize=(8, 12)),
    plt.axes(projection=ccrs.PlateCarree())
)
ax.set_extent((6.4, 6.8, -0.045, .5), crs=ccrs.PlateCarree())

landTuples = (
    ('110m', '#dfe7fdAA', 30), 
    ('50m', '#dfe7fdAA', 30), 
    ('10m', '#3f37c935', 10), 
    ('10m', '#ffffffAA', 2)
)

lands = [
    cfeature.NaturalEarthFeature(
        'physical', 'land', i[0],
        edgecolor=i[1], facecolor='#00000000', linewidth=i[2]
    ) for i in landTuples
]
[ax.add_feature(i, zorder=(-100+(ix+1))) for (ix, i) in enumerate(lands)]

ax.scatter(
    lnd.pointCoords[:, 0], lnd.pointCoords[:, 1],
    color='gray', linestyle='--',
    transform=ccrs.PlateCarree()
)



# import cartopy.io.shapereader as shapereader

# countries = shapereader.natural_earth(resolution='110m',
#                                       category='cultural',
#                                       name='admin_0_countries')


# # Find the China boundary polygon.
# for country in shapereader.Reader(countries).records():
#     if country.attributes['ISO_A3_EH'] == 'MEX':
#         print("Here!")
#         china = country.geometry
