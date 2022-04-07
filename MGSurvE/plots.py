'''Data-Visualization functions.

'''

import matplotlib
from os import path
from math import log
import matplotlib.pyplot as plt
import MGSurvE.constants as cst
import networkx as nx
from sklearn.preprocessing import normalize
import numpy as np


def plotSites(
        fig, ax, 
        sites, pTypes,
        markers=cst.MKRS, colors=cst.MCOL,
        size=350, edgecolors='w', linewidths=1.25,
        zorder=5, transform=None, **kwargs
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
        lineColor='#03045e', 
        lineMin=0.025, lineWidth=20, 
        alphaMin=.5, alphaAmplitude=2.5,
        zorder=0, transform=None, **kwargs
    ):
    """ Plots a transitions matrix.
    
    Parameters:
        fig (matplotlib): Matplotlib fig object.
        ax (matplotlib): Matplotlib ax object.
        transMtx (numpy matrix): Transitions matrix.
        sitesB (numpy array): Coordinates of the vertices origins (sites/traps).
        sitesA (numpy array): Coordinates of the vertices desinations (sites/traps).
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
                lw=max(lineMin, log(1+lineWidth*transMtx[j][i])),
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
        size=350, zorders=(25, -5), fill=True,
        transform=None, 
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
    # if (cNum-tNum) < 0:
    #     raise Exception(
    #         'Less colors ({}) than trap types ({})!'.format(cNum, tNum)
    #     )
    for (i, trap) in enumerate(trapsCoords):
        tType = trapsTypes[i]
        (col, ec) = (colors[tType], edgecolors[0])
        if trapsFixed[i]:
            ec = edgecolors[1]
        transp = 'DD'
        if not fill:
            transp = '00'
        ax.scatter(
            trap[0], trap[1], 
            marker=marker, color=col[:-2]+transp, 
            s=size, zorder=zorders[0],
            edgecolors=ec, linewidths=lws[0]
        )
        for r in trapsKernels[tType]['radii']:
            circle = plt.Circle(
                (trap[0], trap[1]), r, 
                color=col, fill=fill, ls=ls, 
                lw=lws[1], zorder=zorders[1]
            )
            ax.add_patch(circle)
    return (fig, ax)


def plotTrapsNetwork(
        fig, ax,
        transMtx, traps, sites,
        lineColor='#3d0e61', lineWidth=20, 
        alphaMin=.5, alphaAmplitude=2.5,
        zorder=0, transform=None, **kwargs
    ):
    """ Plots the connectivity network of traps in the landscape.

    Parameters:
        fig (matplotlib): Matplotlib fig object.
        ax (matplotlib): Matplotlib ax object.
        transMtx (numpy array): Full transitions matrix.
        traps (numpy array): Traps' coordinates.
        sites (numpy array): Sites' coordinates.
        lineColor (color): Color for the network's line.
        lineWidth (float): Base width for the connections.
        alphaMin (float): Minimum alpha value for connections.
        alphaAmplitude (float): Multiplier for the alpha (proportional to connection).
        zorders (tuple): Z-orders for marker and circles.
        kwargs (dict): Matplotlib's plot-compliant kwargs.
    
    Returns:
        (fig, ax): Matplotlib (fig, ax) tuple.
    """ 
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
        trapsNumber=None, vmin=0, vmax=1, 
        cmap='Purples', linecolor='#222222', linestyle=':', lw=.5,
        ticks=False, 
        **kwargs
    ):
    """ Block matrix plot for the connection network.

    Parameters:
        fig (matplotlib): Matplotlib fig object.
        ax (matplotlib): Matplotlib ax object.
        trapsNumber (int): Number of traps in landscape.
        vmin (float): Lower clipping value.
        vmax (float): Higher clipping value.
        cmap (matplotlib colormap): Matplotlib's colormap object.
        lineColor (color): Color for the block's boundaries.
        lineStyle (matplotlib linestyle): Linestyle for the block matrix' boundaries.
        lw (float): Linewidth for block matrix boundaries.
        ticks (bool): Imshow ticks on/off
        zorders (tuple): Zorders for marker and circles.
        kwargs (dict): Matplotlib's plot-compliant kwargs.
    
    Returns:
        (fig, ax): Matplotlib (fig, ax) tuple.
    """ 
    ax.imshow(
        matrix, 
        cmap=cmap, vmin=vmin, vmax=vmax, aspect='equal',
        **kwargs
    )
    if trapsNumber is not None:
        sitesNumber = matrix.shape[0] - trapsNumber
        plt.axhline(sitesNumber-.5, color=linecolor, ls=linestyle, lw=lw)
        plt.axvline(sitesNumber-.5, color=linecolor, ls=linestyle, lw=lw)
    if not ticks:
        ax.set_xticks([]) 
        ax.set_yticks([]) 
    return (fig, ax)


def plotClean(fig, ax, frame=False, bbox=None, labels=False):
    """ Makes axes equally spaced and removes frame.

    Parameters:
        fig (matplotlib): Matplotlib fig object.
        ax (matplotlib): Matplotlib ax object.
        frame (bool): Flag to remove plot's frame.
    
    Returns:
        (fig, ax): Matplotlib (fig, ax) tuple.
    """ 
    ax.set_aspect('equal')
    if frame is not True:
        ax.axis('off')
    if labels is False:
        ax.set_xticks([])
        ax.set_yticks([])
    if bbox is not None:
        ax.set_xlim(*bbox[0])
        ax.set_ylim(*bbox[1])
    # ax.set_aspect(1/ax.get_data_ratio())
    return (fig, ax)


def plotFitness(
        fig, ax,
        fitness,
        pos=(0.5, 0.5),
        fmt='{:.2f}',
        fontSize=125,
        color='#00000011',
        zorder=5,
        **kwargs
    ):
    """ Adds the fitness value to the plot.

    Parameters:
        fig (matplotlib): Matplotlib fig object.
        ax (matplotlib): Matplotlib ax object.
        pos (floats tuple): Position for the text.
        fmt (string formating): String format for the fitness text.
        fontSize (float): Text's font size.
        color (color): Font color
        zorder (int): Zorder for the text.
        **kwargs: Matplotlib's text kwargs.

    Returns:
        (fig, ax): Matplotlib (fig, ax) tuple.
    """ 
    ax.text(
        pos[0], pos[1], fmt.format(fitness),
        horizontalalignment='center', verticalalignment='center',
        fontsize=fontSize, color=color,
        transform=ax.transAxes, zorder=zorder, **kwargs
    )
    return (fig, ax)


def plotGAEvolution(
        fig, ax,
        gaLog,
        colors={'mean': '#ffffff', 'envelope': '#1565c0'},
        alphas={'mean': .75, 'envelope': 0.5},
        aspect=1/3
    ):
    """ Makes axes equally spaced and removes frame.

    Parameters:
        fig (matplotlib): Matplotlib fig object.
        ax (matplotlib): Matplotlib ax object.
        gaLog (pandas dataframe): Flag to remove plot's frame.
        colors (dict): Mean and envelope colors
        alphas (dict): Mean and envelope alphas
    
    Returns:
        (fig, ax): Matplotlib (fig, ax) tuple.
    """ 
    dta = gaLog
    x = range(dta.shape[0])    
    plt.plot(
        x, dta['avg'], 
        lw=2, color=colors['mean'], alpha=alphas['mean']
    )
    ax.fill_between(
        x, dta['min'], dta['max'], 
        alpha=alphas['envelope'], color=colors['envelope'], lw=0
    )
    ax.set_xlim(0, max(x))
    # ax.set_ylim(0, 5*minFits[-1])
    ax.set_aspect(aspect/ax.get_data_ratio())
    return (fig, ax)


def plotDirectedNetwork(
        fig, ax, 
        sites, pTypes, transMtx,
        markers=cst.MKRS, colors=cst.MCOL,
        alphaNodeMin=1, alphaEdgeMin=1,
        alphaNodeAmplitude=50, alphaEdgeAmplitude=100,
        sizeNodeAmplitude=10000000000, widthEdgeAmplitude=10,
        edgecolors='black', transform=None, **kwargs
    ):
    """ Plots edge and node centrality.
    
    Parameters:
        fig (matplotlib): Matplotlib fig object.
        ax (matplotlib): Matplotlib ax object.
        sites (numpy array): Coordinates of the points.
        pTypes (numpy array): Point types.
        transMtx (numpy matrix): Transitions matrix.
        markers (list): List of marker shapes for point-types (matplotlib).
        colors (list): List of colors for point-types (matplotlib).
        alphaNodeMin (float): Minimum alpha value allowed for nodes.
        alphaEdgeMin (float): Minimum alpha value allowed for edges.
        alphaNodeAmplitude (float): Alpha multiplier for nodes of matrix.
        alphaEdgeAmplitude (float): Alpha multiplier for edges of matrix.
        sizeNodeAmplitude (float): Size multiplier for nodes of matrix.
        widthEdgeAmplitude (float): Width multiplier for edges of matrix.
        edgecolors (color): Edge color.
        kwargs (dict): Matplotlib's plot-compliant kwargs.
    
    Returns:
        (fig, ax): Matplotlib (fig, ax) tuple.
    """
    np.fill_diagonal(transMtx, 0)
    transMtxN = normalize(transMtx, axis=1, norm='l2')
    G = nx.from_numpy_matrix(transMtxN)
    G.remove_edges_from(nx.selfloop_edges(G))

    nodesNum = len(G)
    for i in range(nodesNum):
        keys = G[i]
        for j in range(nodesNum):
            prb = keys.get(j)
            if prb is not None:
                weight = prb['weight']
                if weight > 0:
                    distance = 1 / prb['weight']
                else:
                    distance = np.Inf
                G[i][j]['distance'] = distance

    centrality_nodes = nx.load_centrality(G, weight='distance')
    centrality_edges = nx.edge_betweenness_centrality(G, weight='distance')

    final_G = nx.DiGraph()

    for (i, site) in enumerate(sites):
        final_G.add_node(i, pos=(site[0], site[1]), 
            shape=markers[pTypes[i]], 
            color=colors[pTypes[i]], 
            size=list(centrality_nodes.values())[i])

    for item in centrality_edges.items():
        edge = item[0]
        centrality = item[1]
        final_G.add_edge(int(edge[0]), int(edge[1]), 
            weight=centrality)

    pos = nx.get_node_attributes(final_G, 'pos')

    widths = nx.get_edge_attributes(final_G, 'weight')
    edge_sizes = set(list(widths.values()))

    for shape in set(markers):
        node_list = [node for node in final_G.nodes() if final_G.nodes[node]['shape']==shape]
        nx.draw_networkx_nodes(
            final_G, pos,
            ax=ax,
            nodelist=node_list,
            node_size=[log(1+sizeNodeAmplitude*final_G.nodes[node]['size']) for node in node_list],
            node_color=[final_G.nodes[node]['color'] for node in node_list],
            node_shape=shape,
            alpha=[min(alphaNodeMin, log(1+alphaNodeAmplitude*final_G.nodes[node]['size'])) for node in node_list]
        )

    for width in edge_sizes:
        edge_list = [edge for edge in final_G.edges() if final_G.edges[edge]['weight']==width]
        nx.draw_networkx_edges(
            final_G,pos,
            ax=ax,
            edgelist=edge_list,
            width=log(1+widthEdgeAmplitude*width),
            edge_color=edgecolors,
            alpha=min(alphaEdgeMin, log(1+alphaEdgeAmplitude*width)),
            connectionstyle='arc3,rad=0.2'
        )

    return (fig, ax)


def saveFig(
        fig, ax,
        filepath, filename,
        dpi=300,
        facecolor='w',
        transparent=False,
        bbox_inches='tight',
        pad_inches=0,
        **kwargs
    ):
    """ Save figure to disk.

    Parameters:
        fig (matplotlib): Matplotlib fig object.
        ax (matplotlib): Matplotlib ax object.
        filepath (string): Path for figure export.
        filename (string): Filename for the export.
        dpi (int): Image resolution.
        facecolor (color): Background for the plot.
        transparent (bool): Transparent background.
        bbox_inches (string): Bounding box inches.
        pad_inches (float): Padding inches.
        **kwargs: Matplotlib savefig kwargs.
    
    Returns:
        (fig, ax): Matplotlib (fig, ax) tuple.
    """ 
    fig.savefig(
        path.join(filepath, filename), dpi=dpi,
        facecolor=facecolor, bbox_inches=bbox_inches, 
        pad_inches=pad_inches, transparent=transparent, 
        **kwargs
    )
    return (fig, ax)


def plotTrapsKernels(
        fig, ax, lnd,
        colors=cst.TRP_COLS, maxSca=5, alpha=.75,
        distRange=(0, 100), aspect=.3
    ):
    kers = lnd.trapsKernels
    # dMax = max(max([kers[i]['radii'] for i in range(len(kers))])) * maxSca
    dMax = distRange[1]
    # Generate figure
    for i in range(len(kers)):
        ker = kers[i]
        dists = np.arange(0, dMax+1, dMax/100)
        probs = [ker['kernel'](d, **ker['params']) for d in dists]
        ax.plot(dists, probs, color=colors[i], lw=4, alpha=alpha)
    ax.set_xlim(0, dMax)
    ax.set_ylim(0, 1)
    ax.set_aspect(aspect/ax.get_data_ratio())
    return (fig, ax)


def plotsClearMemory():
    # https://stackoverflow.com/questions/28757348/how-to-clear-memory-completely-of-all-matplotlib-plots
    allfignums = matplotlib.pyplot.get_fignums()
    for i in allfignums:
        fig = matplotlib.pyplot.figure(i)
        fig.clear()
        matplotlib.pyplot.close(fig)


###############################################################################
# Cartopy dependency
###############################################################################
try:
    import cartopy.feature as cfeature
except ImportError:
    print("Cartopy installation was not detected! Geo-boundaries not available!")
else:
    def plotLandBoundary(fig, ax, landTuples=cst.LAND_TUPLES):
        """ Plots the land's boundary as a polygon.

        Parameters:
            fig (matplotlib): Matplotlib fig object.
            ax (matplotlib): Matplotlib ax object.
            landTuples (list of tuples): Check the constants.py file for format.
        
        Returns:
            (fig, ax): Matplotlib (fig, ax) tuple.
        """ 
        lands = [
            cfeature.NaturalEarthFeature(
                'physical', 'land', i[0],
                edgecolor=i[1], facecolor='#00000000', linewidth=i[2]
            ) for i in landTuples
        ]
        [ax.add_feature(i, zorder=(-50+ix)) for (ix, i) in enumerate(lands)]
        return (fig, ax)
