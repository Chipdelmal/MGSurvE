
import math
import time
import random
import numpy as np
import numpy.random as rand
from sklearn.cluster import KMeans
from pointpats import PoissonClusterPointProcess, Window


def ptsRegularGrid(pointsNumber, bbox):
    """ Creates a regular grid of points.
    
    Parameters:
        pointsNumber (int): Number of sites.
        bbox (tuple of tuples): Bounding box in the form ((xLo, xHi), (yLo, yHi)).

    Returns:
        (numpy array): Points' coordinates.
    """
    (xRan, yRan) = bbox
    x = np.linspace(xRan[0], xRan[1], pointsNumber)
    y = np.linspace(yRan[0], yRan[1], pointsNumber)
    coords = np.asarray(np.meshgrid(x, y)).T
    coords = np.concatenate(coords)
    return coords


def ptsDonut(pointsNumber, radii, center=(0, 0)):
    """ Creates a distribution of points laid around a donut shape.
    
    Parameters:
        pointsNumber (int): Number of sites.
        radii (tuple of floats): Minimum and maximum radii for the donut shape.
        center (tuple of floats): Coordinates for the center of the donut (x, y).

    Returns:
        (numpy array): Points' coordinates.
    """
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
    return np.asarray(coords)


def ptsRandUniform(pointsNumber, bbox):
    """ Creates a random unifor distribution of points.
    
    Parameters:
        pointsNumber (int): Number of sites
        bbox (tuple of tuples): Bounding box in the form ((xLo, xHi), (yLo, yHi))

    Returns:
        (numpy array): Points' coordinates
    """
    (xRan, yRan) = bbox
    xy = (
        rand.uniform(*xRan, pointsNumber), 
        rand.uniform(*yRan, pointsNumber)
    )
    coords = list(zip(*xy))
    return np.asarray(coords)


###############################################################################
# Clustering and Aggregating Landscape
###############################################################################

def clusterLandscape(
        pointsCoords, clustersNumber, 
        randomState=time.time(), clusterAlgorithm=KMeans
    ):
    """ .
    
    Parameters:
        pointsCoords (np array): 
        clustersNumber (int):
        randomState (int):
        clusterAlgorithm (function):

    Returns:
        (dict): 
    """
    clObj = clusterAlgorithm(
        n_clusters=clustersNumber,
        random_state=int(randomState)
    )
    clustersObj = clObj.fit(pointsCoords)
    (clusters, centroids) = (
        clustersObj.labels_,
        clustersObj.cluster_centers_
    )
    return {'clusters': clusters, 'centroids': centroids}


def aggregateLandscape(migrationMatrix, clusters):
    """ .
    
    Parameters:
        migrationMatrix (np matrix): 
        clusters (list): 

    Returns:
        (numpy array): 
    """
    matrix_size = len(clusters)
    num_clusters = len(set(clusters))
    aggr_matrix = np.zeros([num_clusters, num_clusters], dtype=float)
    aggr_number = [0]*num_clusters
    for row in range(matrix_size):
        cRow = clusters[row]
        aggr_number[cRow] += 1
        for col in range(matrix_size):
            cCol = clusters[col]
            aggr_matrix[cRow][cCol] += migrationMatrix[row][col]
    for row in range(num_clusters):
        aggr_matrix[row] = [x/aggr_number[row] for x in aggr_matrix[row]]
    return aggr_matrix

def clusterPossion(
    pointsNumber, clustersNumber, 
    radius, randomState=time.time(), 
    bbox=None, polygon=None
    ):
    """ .
    
    Parameters:
        pointsNumber (int): Number of sites
        clustersNumber (int):
        radius (float): Radius of the circle centered on each parent.
        randomState (int):

        bbox (tuple of tuples): Bounding box in the form ((xLo, xHi), (yLo, yHi))
        OR
        polygon (list of tuples): Bounding polygon represented by list of tuples

    Returns:
        (numpy array):  
    """

    if bbox:
        xLo, xHi = bbox[0]
        yLo, yHi = bbox[1]
        parts = [[(xLo, yLo), (xLo, yHi), (xHi, yHi), (xHi, yLo), (xLo, yLo)]]
        window = Window(parts)
    elif polygon:
        window = Window(polygon)

    np.random.seed(randomState)
    csamples = PoissonClusterPointProcess(
        window, pointsNumber, 
        clustersNumber, radius, 1, asPP=True, conditioning=False
    )

    return np.asarray(csamples.realize(pointsNumber))
