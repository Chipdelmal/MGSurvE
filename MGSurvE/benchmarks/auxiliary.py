
import numpy as np
from scipy.interpolate import griddata


def axisRange(x):
    return [min(x), max(x)]


def calcResponseSurface(
    iX, iY, dZ,
    scalers=(1, 1, 1), mthd='linear',
    xAxis='linear', yAxis='linear',
    xLogMin=1e-10, yLogMin=1e-10,
    DXY=(5000, 5000)
):
    (xN, yN, zN) = (
            np.array([float(i/scalers[0]) for i in iX]),
            np.array([float(i/scalers[1]) for i in iY]),
            np.array([float(i/scalers[2]) for i in dZ])
        )
    (xRan, yRan, zRan) = (axisRange(i) for i in (xN, yN, zN))
    # X-Axis scale ------------------------------------------------------------
    if xAxis=='linear':
        xi = np.linspace(xRan[0], xRan[1], DXY[0])
    elif xAxis=='log':
        if xRan[0] == 0:
            xRan[0] = xLogMin
        xi = np.geomspace(xRan[0], xRan[1], DXY[0])
    # Y-Axis scale ------------------------------------------------------------
    if yAxis=='linear':
        yi = np.linspace(yRan[0], yRan[1], DXY[1])
    elif yAxis=='log':
        if yRan[0] == 0:
            yRan[0] =  yLogMin
        yi = np.geomspace(yRan[0], yRan[1], DXY[1])
    # Grid --------------------------------------------------------------------
    zi = griddata((xN, yN), zN, (xi[None, :], yi[:, None]), method=mthd)
    # Return variables
    ranges = (xRan, yRan, zRan)
    grid = (xN, yN, zN)
    surf = (xi, yi, zi)
    return {'ranges': ranges, 'grid': grid, 'surface': surf}


def forceAspect(ax, aspect=1):
    im = ax.get_images()
    extent =  im[0].get_extent()
    ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)