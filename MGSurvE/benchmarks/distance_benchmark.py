import sys
import MGSurvE.auxiliary as aux
import time, random

def main(num = 1000000):
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

if __name__ == '__main__':
    main()
    
