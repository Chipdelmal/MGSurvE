#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import MGSurvE as srv

pts = [
    [0.00, 0.00, 0], 
    [0.25, 0.50, 1], 
    [1.00, 0.15, 0]
]
points = pd.DataFrame(pts, columns=['x', 'y', 't'])

msk = [
    [.6, .4],
    [.3, .7]
]

lnd = srv.Landscape(points, pointTypesMask=msk)
# lnd.calculatePointsDistances()
# lnd.distanceMatrix
lnd.distanceMatrix
lnd.migrationMatrix
lnd.maskedMigration

