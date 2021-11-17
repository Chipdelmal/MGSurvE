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

MCOL = ('#e0c3fc', '#bdb2ff', '#a0c4ff', '#ffd6a5', '#caffbf', '#d0d1ff')
"""A cute pastel colors list."""

PINK_NAVY = col.colorPaletteFromHexList(['#e0c3fc',  '#00296b'])
"""Pink to Navy Blue cmap."""

TRP_COLS = {
    0: '#f7258511', 1: '#9381ff11', 2: '#5ddeb111', 
    3: '#f038ff11', 4: '#fe5f5511', 5: '#e2ef7011', 
}
"""Base colors for trap types."""