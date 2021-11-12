
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
        plt.scatter(
            site[0], site[1], 
            marker=markers[pTypes[i]], color=colors[pTypes[i]], 
            s=size, zorder=zorder, edgecolors=edgecolors, linewidths=linewidths
        )
    return (fig, ax)


def plotNetwork(
        fig, ax, 
        transMtx, sites, 
        lineColor='#03045e', lineWidth=20, 
        alphaMin=.5, alphaAmplitude=2.5,
        zorder=0, **kwargs
    ):
    """ Plots a transitions matrix.
    
    Parameters:
        fig (matplotlib): Matplotlib fig object.
        ax (matplotlib): Matplotlib ax object.
        transMtx (numpy matrix): Transitions matrix.
        sites (numpy array): Coordinates of the vertices.
        lineColor (color): Color for the network.
        lineWidth (float): Amplitude for the linewidth.
        alphaMin (float): Minimum alpha value allowed.
        alphaAmplitude (float): Alpha multiplier for matrix.
        zorder (int): Matplotlib's zorder.
        kwargs (dict): Matplotlib's plot-compliant kwargs.
    
    Returns:
        (fig, ax): Matplotlib (fig, ax) tuple.
    """
    (aNum, bNum) = (sites.shape[0], sites.shape[0])
    for j in range(aNum):
        src = sites[j]
        for i in range(bNum):
            snk = sites[i]
            plt.plot(
                [src[0], snk[0]], [src[1], snk[1]],
                lw=log(1+lineWidth*transMtx[j][i]),
                alpha=min(alphaMin, log(1+alphaAmplitude*transMtx[j][i])),
                c=lineColor, zorder=zorder,
                **kwargs
            )
    return (fig, ax)


def plotTraps(
        fig, ax,
        trapsCoords, trapsTypes, radii,
        colors=['#f72585FA', ], markers=["X", ],
        **kwargs
    ):
    """ Plots the traps with the radii of effectiveness.

    """
    # for trap in trapsCoords:
    #     plt.scatter(
    #         trap[0], trap[1], 
    #         marker="X", color=color, s=600, zorder=25,
    #         edgecolors='w', linewidths=2
    #     )
    #     for r in radii:
    #         circle = plt.Circle(
    #             (trap[0], trap[1]), r, 
    #             color='#f7258509', fill=True, ls=':', lw=0, zorder=0
    #         )
    #         ax.add_patch(circle)
    return (fig, ax)