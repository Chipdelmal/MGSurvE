import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import MGSurvE.auxiliary as aux
import time, random

def benchmark(num = 1000000):
    #Set random seed for consistent benchmarks
    random.seed (0x7126434a2ea2a259e9f4196cbb343b1e6d4c2fc8)

###############################################################################
# Create pairs of random points
###############################################################################
    # Longitude then latitude
    def random_point():
        return [random.uniform(-180, 180), random.uniform(-90, 90)]

    def random_index():
        return random.randint(0, num - 1)

    points = [random_point() for i in range(num)]
    rand_points = [points[random_index()] for i in range(num)]
###############################################################################
# Test Vincenty
###############################################################################
    vin_beg_time = time.perf_counter()
    vin_dist = 0
    for (pt0, pt1) in zip(points, rand_points):
        vin_dist += aux.vincentyDistance(pt0, pt1)
    vin_end_time = time.perf_counter()

    vin_elap_time = (vin_end_time - vin_beg_time)

    print("Vincenity Performance:")
    print(f"Completed in {vin_elap_time} seconds")
    print("\n")
###############################################################################
# Test Cheap Ruler
###############################################################################
    cheap_beg_time = time.perf_counter()
    cheap_dist = 0
    for (pt0, pt1) in zip(points, rand_points):
        cheap_dist += aux.cheapRuler(pt0, pt1)
    cheap_end_time = time.perf_counter()

    cheap_elap_time = (cheap_end_time - cheap_beg_time)

    print("Cheap Ruler Performance:")
    print(f"Completed in {cheap_elap_time} seconds")
    print("\n")
###############################################################################
# Test Haversine
###############################################################################
    haver_beg_time = time.perf_counter()
    haver_dist = 0
    for (pt0, pt1) in zip(points, rand_points):
        haver_dist += aux.haversineDistance(pt0, pt1)
    haver_end_time = time.perf_counter()

    haver_elap_time = (haver_end_time - haver_beg_time)

    print("Haversine Performance:")
    print(f"Completed in {haver_elap_time} seconds")
    print("\n")

def dataframes(num = 1000000):
    #Set random seed for consistent benchmarks
    random.seed (0x7126434a2ea2a259e9f4196cbb343b1e6d4c2fc8)

    fig, axes = plt.subplots(2, 3)
###############################################################################
# Create pairs of random points
###############################################################################
    # Longitude then latitude
    def random_point():
        return [random.uniform(-180, 180), random.uniform(-90, 90)]

    def random_index():
        return random.randint(0, num - 1)

    points = [random_point() for i in range(num)]
    rand_points = [points[random_index()] for i in range(num)]
###############################################################################
# Test Vincenty
###############################################################################
    pt0Lon = []
    pt0Lat = []
    pt1Lon = []
    pt1Lat = []
    distances = []
    times = []

    for (pt0, pt1) in zip(points, rand_points):
        beg_time = time.perf_counter()
        dist = aux.vincentyDistance(pt0, pt1)
        end_time = time.perf_counter()

        elap_time = (end_time - beg_time)

        pt0Lon.append(pt0[0])
        pt0Lat.append(pt0[1])
        pt1Lon.append(pt1[0])
        pt1Lat.append(pt1[1])
        distances.append(dist)
        times.append(elap_time)

    d = {'pt0 Longitude' : pt0Lon, 'pt0 Latitude' : pt0Lat,
    'pt1 Longitude' : pt1Lon, 'pt1 Latitude' : pt1Lat,
    'Calculated distance' : distances, 'Time (seconds)' : times}

    vinDF = pd.DataFrame(data = d)

###############################################################################
# Test Cheap Ruler
###############################################################################
    pt0Lon = []
    pt0Lat = []
    pt1Lon = []
    pt1Lat = []
    distances = []
    times = []

    for (pt0, pt1) in zip(points, rand_points):
        beg_time = time.perf_counter()
        dist = aux.cheapRuler(pt0, pt1)
        end_time = time.perf_counter()

        elap_time = (end_time - beg_time)

        pt0Lon.append(pt0[0])
        pt0Lat.append(pt0[1])
        pt1Lon.append(pt1[0])
        pt1Lat.append(pt1[1])
        distances.append(dist)
        times.append(elap_time)

    d = {'pt0 Longitude' : pt0Lon, 'pt0 Latitude' : pt0Lat,
    'pt1 Longitude' : pt1Lon, 'pt1 Latitude' : pt1Lat,
    'Calculated distance' : distances, 'Time (seconds)' : times}

    cheapDF = pd.DataFrame(data = d)
###############################################################################
# Test Haversine
###############################################################################
    pt0Lon = []
    pt0Lat = []
    pt1Lon = []
    pt1Lat = []
    distances = []
    times = []

    for (pt0, pt1) in zip(points, rand_points):
        beg_time = time.perf_counter()
        dist = aux.haversineDistance(pt0, pt1)
        end_time = time.perf_counter()

        elap_time = (end_time - beg_time)

        pt0Lon.append(pt0[0])
        pt0Lat.append(pt0[1])
        pt1Lon.append(pt1[0])
        pt1Lat.append(pt1[1])
        distances.append(dist)
        times.append(elap_time)

    d = {'pt0 Longitude' : pt0Lon, 'pt0 Latitude' : pt0Lat,
    'pt1 Longitude' : pt1Lon, 'pt1 Latitude' : pt1Lat,
    'Calculated distance' : distances, 'Time (seconds)' : times}

    havDF = pd.DataFrame(data = d)
###############################################################################
# Save dataframes
###############################################################################
    vinDF.to_csv('vincenty_distances.csv', index=False)
    cheapDF.to_csv('cheapRuler_distances.csv', index=False)
    havDF.to_csv('haversine_distances.csv', index=False)

def violinPlots(num = 1000000):
    #Set random seed for consistent benchmarks
    random.seed (0x7126434a2ea2a259e9f4196cbb343b1e6d4c2fc8)
    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(9, 4), sharey=True)

###############################################################################
# Helper functions
###############################################################################
    def is_outlier(points, thresh=3.5):
        if len(points.shape) == 1:
            points = points[:,None]
        median = np.median(points, axis=0)
        diff = np.sum((points - median)**2, axis=-1)
        diff = np.sqrt(diff)
        med_abs_deviation = np.median(diff)

        modified_z_score = 0.6745 * diff / med_abs_deviation

        return modified_z_score > thresh

###############################################################################
# Create pairs of random points
###############################################################################
    # Longitude then latitude
    def random_point():
        return [random.uniform(-180, 180), random.uniform(-90, 90)]

    def random_index():
        return random.randint(0, num - 1)

    points = [random_point() for i in range(num)]
    rand_points = [points[random_index()] for i in range(num)]
###############################################################################
# Test Vincenty
###############################################################################
    times = []

    for (pt0, pt1) in zip(points, rand_points):
        beg_time = time.perf_counter()
        dist = aux.vincentyDistance(pt0, pt1)
        end_time = time.perf_counter()

        elap_time = (end_time - beg_time)
        times.append(elap_time*1000)

    times = np.array(times)
    times = times[~is_outlier(times)]   

    axes[0].set_title("Vincenty Execution Times")
    axes[0].set_ylabel("Times (milliseconds)")
    axes[0].violinplot(times)
    axes[0].xaxis.set_tick_params(direction='out')
    axes[0].xaxis.set_ticks_position('bottom')
    axes[0].set_xticks(np.arange(1, 2), labels=['Vincenty'])
    axes[0].set_xlim(0.25, 1.75)

###############################################################################
# Test Cheap Ruler
###############################################################################
    times = []

    for (pt0, pt1) in zip(points, rand_points):
        beg_time = time.perf_counter()
        dist = aux.cheapRuler(pt0, pt1)
        end_time = time.perf_counter()

        elap_time = (end_time - beg_time)
        times.append(elap_time*1000)

    times = np.array(times) 
    times = times[~is_outlier(times)]     

    axes[1].set_title("Cheap Ruler Execution Times")
    axes[1].violinplot(times)
    axes[1].xaxis.set_tick_params(direction='out')
    axes[1].xaxis.set_ticks_position('bottom')
    axes[1].set_xticks(np.arange(1, 2), labels=['Cheap Ruler'])
    axes[1].set_xlim(0.25, 1.75)
###############################################################################
# Test Haversine
###############################################################################
    times = []

    for (pt0, pt1) in zip(points, rand_points):
        beg_time = time.perf_counter()
        dist = aux.haversineDistance(pt0, pt1)
        end_time = time.perf_counter()

        elap_time = (end_time - beg_time)
        times.append(elap_time*1000)

    times = np.array(times)   
    times = times[~is_outlier(times)]   

    axes[2].set_title("Haversine Execution Times")
    axes[2].violinplot(times)
    axes[2].xaxis.set_tick_params(direction='out')
    axes[2].xaxis.set_ticks_position('bottom')
    axes[2].set_xticks(np.arange(1, 2), labels=['Haversine'])
    axes[2].set_xlim(0.25, 1.75)
###############################################################################
# Save Figure
###############################################################################
    plt.savefig("distancePlots5.png")

def errorPlots(num = 100):
    #Set random seed for consistent benchmarks
    random.seed (0x7126434a2ea2a259e9f4196cbb343b1e6d4c2fc8)
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(9, 4), sharey=True)

###############################################################################
# Create pairs of random points
###############################################################################
    # Longitude then latitude
    def random_point():
        return [random.uniform(-180, 180), random.uniform(-90, 90)]

    def random_index():
        return random.randint(0, num - 1)

    points = [random_point() for i in range(num)]
    rand_points = [points[random_index()] for i in range(num)]
###############################################################################
# Test Vincenty
###############################################################################
    vincentyDistances = []

    for (pt0, pt1) in zip(points, rand_points):
        dist = aux.vincentyDistance(pt0, pt1)

        vincentyDistances.append(dist)

    vincentyDistances = np.array(vincentyDistances)

###############################################################################
# Test Cheap Ruler
###############################################################################
    cheapDistances = []

    for (pt0, pt1) in zip(points, rand_points):
        dist = aux.cheapRuler(pt0, pt1)

        cheapDistances.append(dist)

    cheapDistances = np.array(cheapDistances)

    axes[0].set_title("Vincenty vs Cheap Ruler")
    axes[0].scatter(x=vincentyDistances, y=[cheapDistances[i] - vincentyDistances[i] for i in range(len(vincentyDistances)) if vincentyDistances[i] != 0])
    axes[0].set_xlabel("Vincenty Distance (meters)")
    axes[0].set_ylabel("(Cheap Ruler - Vincenty) Distance Difference")
###############################################################################
# Test Haversine
###############################################################################
    haversineDistances = []

    for (pt0, pt1) in zip(points, rand_points):
        dist = aux.haversineDistance(pt0, pt1)

        haversineDistances.append(dist)

    haversineDistances = np.array(haversineDistances)

    axes[1].set_title("Vincenty vs Haversine")
    axes[1].scatter(x=vincentyDistances, y=[haversineDistances[i] - vincentyDistances[i] for i in range(len(vincentyDistances)) if vincentyDistances[i] != 0])
    axes[1].set_xlabel("Vincenty Distance (meters)")
    axes[1].set_ylabel("(Haversine - Vincenty) Distance Difference")
###############################################################################
# Save Figure
###############################################################################
    plt.savefig("errorPlots.png")

if __name__ == '__main__':
    violinPlots()
    
