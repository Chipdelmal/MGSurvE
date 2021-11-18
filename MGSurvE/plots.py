
from math import log
import matplotlib.pyplot as plt
import MGSurvE.constants as cst


def plotSites(
        fig, ax, 
        sites, pTypes,
        markers=cst.MKRS, colors=cst.MCOL,
        size=150, edgecolors='w', linewidths=1.25,
        zorder=5, **kwargs
    ):
    """ Plots a transitions matrix.
    
    Parameters:
        fig (matplotlib): Matplotlib fig object.
        ax (matplotlib): Matplotlib ax object.
        sites (numpy array): Coordinates of the points.
        pTypes (numpy array): Point types.
        markers (list): List of marker shapes for point-types (matplotlib).
        colors (list): List of colors for point-types (matplotlib).
        size (float): Marker size.
        edgecolors (color): Edge color for markers.
        linewidths (float): Edge line width for markers.
        zorder (int): Matplotlib's zorder.
        kwargs (dict): Matplotlib's plot-compliant kwargs.
    
    Returns:
        (fig, ax): Matplotlib (fig, ax) tuple.
    """
    for (i, site) in enumerate(sites):
        ax.scatter(
            site[0], site[1], 
            marker=markers[pTypes[i]], color=colors[pTypes[i]], 
            s=size, zorder=zorder, edgecolors=edgecolors, linewidths=linewidths
        )
    return (fig, ax)


def plotMigrationNetwork(
        fig, ax, 
        transMtx, sitesB, sitesA,
        lineColor='#03045e', lineWidth=20, 
        alphaMin=.5, alphaAmplitude=2.5,
        zorder=0, **kwargs
    ):
    """ Plots a transitions matrix.
    
    Parameters:
        fig (matplotlib): Matplotlib fig object.
        ax (matplotlib): Matplotlib ax object.
        transMtx (numpy matrix): Transitions matrix.
        sitesB (numpy array): Coordinates of the vertices origins (sites/traps).
        sitesA (numpy array): Coordinates of the vertices origins (sites).
        lineColor (color): Color for the network.
        lineWidth (float): Amplitude for the linewidth.
        alphaMin (float): Minimum alpha value allowed.
        alphaAmplitude (float): Alpha multiplier for matrix.
        zorder (int): Matplotlib's zorder.
        kwargs (dict): Matplotlib's plot-compliant kwargs.
    
    Returns:
        (fig, ax): Matplotlib (fig, ax) tuple.
    """
    (aNum, bNum) = (sitesA.shape[0], sitesB.shape[0])
    for j in range(aNum):
        src = sitesA[j]
        for i in range(bNum):
            snk = sitesB[i]
            ax.plot(
                [src[0], snk[0]], [src[1], snk[1]],
                lw=log(1+lineWidth*transMtx[j][i]),
                alpha=min(alphaMin, log(1+alphaAmplitude*transMtx[j][i])),
                c=lineColor, zorder=zorder,
                **kwargs
            )
    return (fig, ax)


def plotTraps(
        fig, ax,
        trapsCoords, trapsTypes, trapsKernels, trapsFixed,
        colors=cst.TRP_COLS, marker="X",
        edgecolors=('w', 'k'), lws=(2, 0), ls=':',
        size=200, zorders=(25, -5),
        **kwargs
    ):
    """ Plots the traps with the radii of effectiveness.

    Parameters:
        fig (matplotlib): Matplotlib fig object.
        ax (matplotlib): Matplotlib ax object.
        trapsCoords (numpy array): Coordinates of the vertices.
        trapsTypes (list ints): Trap types IDs.
        trapsKernels (dict): Dictionary of traps kernels.
        colors (dict): List of colors for different trap types.
        marker (mrk): Marker type for matplotlib.
        edgecolor (color): Edgecolor for trap marker.
        lws (tuple): Line widths for marker and radii (in order).
        ls (str): Linestyle for the radii.
        size (float): Size of the marker.
        zorders (tuple): Zorders for marker and circles.
        kwargs (dict): Matplotlib's plot-compliant kwargs.
    
    Returns:
        (fig, ax): Matplotlib (fig, ax) tuple.

    """
    (cNum, tNum) = (len(colors), len(trapsTypes))
    if (cNum-tNum) < 0:
        raise Exception(
            'Less colors ({}) than trap types ({})!'.format(cNum, tNum)
        )
    for (i, trap) in enumerate(trapsCoords):
        tType = trapsTypes[i]
        (col, ec) = (colors[tType], edgecolors[0])
        if trapsFixed[i]:
            ec = edgecolors[1]
        ax.scatter(
            trap[0], trap[1], 
            marker=marker, color=col[:-2]+'DD', 
            s=size, zorder=zorders[0],
            edgecolors=ec, linewidths=lws[0]
        )
        for r in trapsKernels[tType]['radii']:
            circle = plt.Circle(
                (trap[0], trap[1]), r, 
                color=col, fill=True, ls=ls, 
                lw=lws[1], zorder=zorders[1]
            )
            ax.add_patch(circle)
    return (fig, ax)


def plotTrapsNetwork(
        fig, ax,
        transMtx, traps, sites,
        lineColor='#3d0e61', lineWidth=20, 
        alphaMin=.5, alphaAmplitude=2.5,
        zorder=0, **kwargs
    ):
    ptsNum = sites.shape[0]
    txMtx = transMtx[:ptsNum, ptsNum:]
    (fig, ax) = plotMigrationNetwork(
        fig, ax, 
        txMtx, traps, sites,
        lineColor=lineColor, lineWidth=lineWidth,
        alphaMin=alphaMin, alphaAmplitude=alphaAmplitude,
        zorder=zorder,
        **kwargs
    )
    return (fig, ax)


def plotMatrix(
        fig, ax,
        matrix, 
        trapsNumber=None, vmin=0, vmax=None, 
        cmap='Purples', linecolor='#222222', linestyle=':', ticks=False, lw=.25,
        **kwargs
    ):
    ax.imshow(
        matrix, 
        cmap=cmap, vmin=vmin, vmax=vmax, aspect='equal',
        **kwargs
    )
    sitesNumber = matrix.shape[0] - trapsNumber
    if trapsNumber is not None:
        plt.axhline(sitesNumber-.5, color=linecolor, ls=linestyle, lw=lw)
        plt.axvline(sitesNumber-.5, color=linecolor, ls=linestyle, lw=lw)
    if not ticks:
        ax.set_xticks([]) 
        ax.set_yticks([]) 
    return (fig, ax)


def plotClean(fig, ax, frame=False):
    ax.set_aspect('equal')
    if frame is not True:
        ax.axis('off')
    return (fig, ax)