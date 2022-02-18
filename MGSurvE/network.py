'''Network-analysis operations.

'''

import numpy as np
import networkx as nx
from sklearn.preprocessing import normalize


def calculateNetworkCentralities(
        transitionsMatrix,
        nodeCentrality=nx.load_centrality,
        edgeCentrality=nx.edge_betweenness_centrality,
        weight='distance'
    ):
    # Setup network object ----------------------------------------------------
    transMtx = np.copy(transitionsMatrix)
    np.fill_diagonal(transMtx, 0)
    transMtxN = normalize(transMtx, axis=1, norm='l2')
    G = nx.from_numpy_matrix(transMtxN)
    G.remove_edges_from(nx.selfloop_edges(G))
    # Add distances to network object -----------------------------------------
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
    # Get centralities --------------------------------------------------------
    centrality_nodes = nodeCentrality(G, weight=weight)
    centrality_edges = edgeCentrality(G, weight=weight)
    return {'edges': centrality_edges, 'nodes': centrality_nodes, 'graph': G}
