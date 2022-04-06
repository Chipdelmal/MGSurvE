
import itertools
import numpy as np
import pandas as pd
from copy import deepcopy
import matplotlib.pyplot as plt
import MGSurvE as srv

LND = 2
ptsNum = 20
bbox = ((-10, 10), (-10, 10))
radii = (20, 25)
###############################################################################
# Select landscape type
###############################################################################
if LND == 0:
    xy = srv.ptsRegularGrid(ptsNum, bbox).T
elif LND == 1:
    xy = srv.ptsDonut(ptsNum, radii).T
elif LND == 2:
    xy = srv.ptsRandUniform(ptsNum, bbox).T
###############################################################################
# Generate landscape object
###############################################################################
points = pd.DataFrame({'x': xy[0], 'y': xy[1], 't': [0]*xy.shape[1]})
lnd = srv.Landscape(
    points, attractionVector=[1]*ptsNum
)
lndOrg = deepcopy(lnd)
orgMig = lndOrg.maskedMigration.T
###############################################################################
# Plot
###############################################################################
(fig, ax) = plt.subplots(1, 1, figsize=(15, 15), sharey=False)
(fig, ax) = lnd.plotSites(fig, ax)
(fig, ax) = lnd.plotMaskedMigrationNetwork(fig, ax)
for (i, coord) in enumerate(lnd.pointCoords):
    ax.text(coord[0], coord[1], str(i), zorder=100, ha='center', va='center')






(HI, LO) = (100, 5)
lnd = srv.Landscape(
    points, 
    attractionVector=([HI]*(ptsNum//2)) +([LO]*(ptsNum//2))
    # list(itertools.chain(*[[HI]*int((ptsNum**2/2))+([LO]*int(ptsNum**2/2))]))
)
modMig = lnd.maskedMigration.T
###############################################################################
# Plot
###############################################################################
(fig, ax) = plt.subplots(1, 1, figsize=(15, 15), sharey=False)
(fig, ax) = lnd.plotSites(fig, ax)
(fig, ax) = lnd.plotMaskedMigrationNetwork(fig, ax)
for (i, coord) in enumerate(lnd.pointCoords):
    ax.text(coord[0], coord[1], str(i), zorder=100, ha='center', va='center')




difference = not(all(np.isclose(np.sum(orgMig, axis=1), np.sum(modMig, axis=1))))
consistency = all(np.isclose(np.sum(orgMig, axis=0), np.sum(modMig, axis=0)))
(difference, consistency)

# In-degree of attractive nodes (HI)
(bse, tst) = (
    np.sum(orgMig[:ptsNum//2, :ptsNum//2], axis=1),
    np.sum(modMig[:ptsNum//2, :ptsNum//2], axis=1)
)
iHi = sum([b < t for (b, t) in zip(bse, tst)])
# In-degree of non-attractive nodes (LO)
(bse, tst) = (
    np.sum(orgMig[ptsNum//2:, ptsNum//2:], axis=1),
    np.sum(modMig[ptsNum//2:, ptsNum//2:], axis=1)
)
iLo = sum([b > t for (b, t) in zip(bse, tst)])

# Out-degree of attractive nodes (HI)
bse = np.sum(orgMig[:ptsNum//2, ptsNum//2:], axis=0)
tst = np.sum(modMig[:ptsNum//2, ptsNum//2:], axis=0)
sum([t > b for (b, t) in zip(bse, tst)])

[t < b for (b, t) in zip(bse, tst)]

# Out-degree of attractive nodes (LO)
bse = np.sum(orgMig[ptsNum//2:, ptsNum//2:], axis=0)
tst = np.sum(modMig[ptsNum//2:, ptsNum//2:], axis=0)
oLo = sum([t < b for (b, t) in zip(bse, tst)])



