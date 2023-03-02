
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sympy.solvers import solve, nsolve
from sympy import Symbol
import MGSurvE as srv

LND = 2
ptsNum = 100
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
    

radius = 200    

pointsNumber = 100
step = 50



(x, y, r, nodesN) = ([0], [0], step, 4)
for rad in np.arange(step, radius+1, step):
    # Calculate angles for nodes
    angleDelta = (2*math.pi)/nodesN
    angles = np.arange(0, 2*math.pi, angleDelta)
    # calculating coordinates
    (xN, yN) = (
        [rad*math.cos(a) for a in angles],
        [rad*math.sin(a) for a in angles]
    )
    x.extend(xN)
    y.extend(yN)
    print(f'{len(angles)} - {len(xN)} - {nodesN} - {len(x)}')
    # Update variables
    r = r+step
    nodesN = int(nodesN*2)
print(len(x))
# print(nodesN/2+5)
delta = radius/step
print(2*(2**(delta+1))-3)

points = pd.DataFrame({'x': x, 'y': y, 't': [0]*len(x)})
lnd = srv.Landscape(points)

(fig, ax) = plt.subplots(1, 1, figsize=(15, 15), sharey=False)
lnd.plotSites(fig, ax)
# lnd.plotMigrationNetwork(fig, ax)
srv.plotClean(fig, ax, frame=False)
###############################################################################
# Solving equation
###############################################################################
p = len(x)
r = radius
(r*math.log(2))/(math.log((p+3)/8))



r = 200    
pN = 200

s = Symbol('s')
step = float(nsolve(2*(2**(r/s+1))-3-pN, s, 10))

delta = radius/step

(x, y, r, nodesN) = ([0], [0], step, 4)
for rad in np.arange(step, radius+1, step):
    # Calculate angles for nodes
    angleDelta = (2*math.pi)/nodesN
    angles = np.arange(0, 2*math.pi, angleDelta)
    # calculating coordinates
    (xN, yN) = (
        [rad*math.cos(a) for a in angles],
        [rad*math.sin(a) for a in angles]
    )
    x.extend(xN)
    y.extend(yN)
    print(f'{len(angles)} - {len(xN)} - {nodesN} - {len(x)}')
    # Update variables
    r = r+step
    nodesN = int(nodesN*2)
    
len(x)


def ptsRegularCircle(pointsNumber, radius, solStart=10):
    (pN, s, r) = (pointsNumber, Symbol('s'), radius)
    # Solve for radius step size for approximate number of points
    step = int(nsolve(2*(2**(r/s+1))-3-pN, s, solStart))
    # Iterate through generated points (origin excluded)
    (x, y, r, nodesN) = ([0], [0], step, 4)
    for rad in np.arange(step, radius+1, step):
        # Calculate angles for nodes
        angleDelta = (2*math.pi)/nodesN
        angles = np.arange(0, 2*math.pi, angleDelta)
        # Map to xy coordinates
        (xN, yN) = (
            [rad*math.cos(a) for a in angles],
            [rad*math.sin(a) for a in angles]
        )
        x.extend(xN)
        y.extend(yN)
        # Update variables
        r = r+step
        nodesN = int(nodesN*2)
    # Convert to array and return
    xy = np.array([x, y])
    return xy

xy = ptsRegularCircle(100, 200)
points = pd.DataFrame({'x': xy[0], 'y': xy[1], 't': [0]*xy.shape[1]})
lnd = srv.Landscape(points)

(fig, ax) = plt.subplots(1, 1, figsize=(15, 15), sharey=False)
lnd.plotSites(fig, ax)
###############################################################################
# Generate landscape object
###############################################################################
points = pd.DataFrame({'x': xy[0], 'y': xy[1], 't': [0]*xy.shape[1]})
lnd = srv.Landscape(points)
###############################################################################
# Plot
###############################################################################
(fig, ax) = plt.subplots(1, 1, figsize=(15, 15), sharey=False)
lnd.plotSites(fig, ax)
lnd.plotMigrationNetwork(fig, ax)
srv.plotClean(fig, ax, frame=False)