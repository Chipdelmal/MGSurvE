'''Constants used across the pkg (colors, symbols, bio).

'''

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math
import MGSurvE.colors as col

###############################################################################
# Bio
###############################################################################
AEDES_EXP_PARAMS = [0.01848777, 1.0e-10, math.inf]
"""Aedes aegypti's migration parameters for kernel."""

SHORT_EXP_PARAMS = [1, 1.0e-10, math.inf]
"""Dummy migration parameters for short-distance kernel."""

MEDIUM_MOV_EXP_PARAMS = [.075, 1.0e-10, math.inf]
"""Dummy migration parameters for short-distance kernel."""

###############################################################################
# Traps
###############################################################################
BASIC_EXP_TRAP = {'A': 0.5, 'b': 0.15}
"""Generic params for an exponential-decay trap"""

###############################################################################
# Plots
###############################################################################
MKRS = ('o', '^', 's', 'p', 'd', 'X')
"""Markers for point-types"""

MCOL = ('#bdb2ff', '#a0c4ff', '#e0c3fc', '#ffd6a5', '#caffbf', '#d0d1ff')
"""A cute pastel colors list."""

PINK_NAVY = col.colorPaletteFromHexList(['#e0c3fc',  '#00296b'])
"""Pink to Navy Blue cmap."""

TRP_COLS = {
    0: '#f7258515', 1: '#5ddeb125', 2: '#fe5f5515', 
    3: '#f038ff15', 4: '#e2ef7015', 5: '#9381ff15', 
}
"""Base colors for trap types."""

LAND_TUPLES = (
    # ('110m', '#dfe7fdAA', 30), ('50m', '#dfe7fdAA', 30), 
    ('10m', '#dfe7fd55', 30), ('10m', '#ffffffDD', 5)
)
"""Base colors for land boundaries."""