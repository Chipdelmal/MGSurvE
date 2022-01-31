from pointpats import PoissonClusterPointProcess, Window, poly_from_bbox, PointPattern
import matplotlib.pyplot as plt
import numpy as np

parts = [[(0.0, 0.0), (0.0, 10.0), (10.0, 10.0), (10.0, 0.0), (0.0, 0.0)]]
window = Window(parts)

np.random.seed(5)
csamples = PoissonClusterPointProcess(window, 200, 10, 0.5, 1, asPP=True, conditioning=False)

pp_pcp = csamples.realizations[0]
pp_pcp.plot(window=True, hull=True, title='Rectangle Clustered Point Pattern')

plt.savefig('test_poly.jpeg')