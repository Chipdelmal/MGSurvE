from pointpats import PoissonClusterPointProcess, Window, poly_from_bbox, PointPattern
import matplotlib.pyplot as plt
import numpy as np
import MGSurvE as srv

parts = [[(0.0, 0.0), (0.0, 10.0), (10.0, 10.0), (10.0, 0.0), (0.0, 0.0)]]
window = Window(parts)

np.random.seed(5)
csamples = PoissonClusterPointProcess(window, 200, 10, 0.5, 1, asPP=True, conditioning=False)

pp_pcp = csamples.realizations[0]
pp_pcp.plot(window=True, hull=True, title='Rectangle Clustered Point Pattern')

plt.savefig('test_poly.png')

tst = srv.clusterPossion(
    100, 5, 100, 
    polygon="/home/chipdelmal/Documents/WorkSims/MGSurvE_Benchmarks/Comoros/GEO/com_admbnda_adm1_cosep_ocha_20191205.shp"
)