#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math
import MGSurvE.colors as col

###############################################################################
# Bio
###############################################################################
AEDES_EXP_PARAMS = [0.01848777, 1.0e-10, math.inf]
"""Aedes aegypti's migration parameters for kernel."""


###############################################################################
# Plots
###############################################################################
MKRS = ('o', '^', 's', 'p', 'd', 'X')
"""Markers for point-types"""

MCOL = ('#e0c3fc', '#bdb2ff', '#a0c4ff', '#ffd6a5', '#caffbf')
"""A cute pastel colors list."""

PINK_NAVY = col.colorPaletteFromHexList(['#e0c3fc',  '#00296b'])
"""Pink to Navy Blue cmap."""

