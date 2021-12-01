
import math
import random
import numpy as np
import numpy.random as rand


def ptsRegularGrid(pointsNumber, bbox):
    (xRan, yRan) = bbox
    x = np.linspace(xRan[0], xRan[1], pointsNumber)
    y = np.linspace(yRan[0], yRan[1], pointsNumber)
    coords = np.asarray(np.meshgrid(x, y)).T
    coords = np.concatenate(coords)
    return coords


def ptsDonut(pointsNumber, radii, center):
    coords = []
    for i in range(pointsNumber):
        radius = (radii[0], radii[1])
        angle = 2 * math.pi * random.random()
        # random radius
        r = rand.uniform(radius[0], radius[1])
        # calculating coordinates
        x = r * math.cos(angle) + center[0]
        y = r * math.sin(angle) + center[1]
        coords.append([x, y])
    return coords


def ptsRandUniform(pointsNumber, bbox):
    (xRan, yRan) = bbox
    xy = (
        rand.uniform(*xRan, pointsNumber), 
        rand.uniform(*yRan, pointsNumber)
    )
    coords = list(zip(xy))
    return coords