import sys
import MGSurvE.auxiliary as aux
import time, random

def main(num = 1000000):
    random.seed (0x7126434a2ea2a259e9f4196cbb343b1e6d4c2fc8)

    # Longitude then latitude
    def random_point():
        return (random.uniform(-180, 180), random.uniform(-90, 90))

    def random_index():
        return random.randint(0, num - 1)

    points = [random_point for i in range(num)]
    rand_points = [points[random_index()] for i in range(num)]

    vin_beg_time = time.perf_counter()
    vin_dist = 0
    for (pt0, pt1) in zip(points, rand_points):
        vin_dist += aux.vincentyDistance(pt0, pt1)
    vin_end_time = time.perf_counter()

    vin_elap_time = (vin_end_time - vin_beg_time)

    print("Vincinity Performance:")
    print(vin_elap_time)
    print("\n")

if __name__ == '__main__':
    main()
    
