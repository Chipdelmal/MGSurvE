#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import MGSurvE as srv

xy = np.asarray([[0, 0], [1, 0]])
pTypes = [0, 1]
pTypesMask = np.asarray([[0]])

lnd = srv.Landscape(xy, pointTypes=pTypes, pointTypesMask=pTypesMask)
# lnd.calculatePointsDistances()
lnd.distanceMatrix