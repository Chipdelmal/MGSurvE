#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap

def colorPaletteFromHexList(clist):
    c = mcolors.ColorConverter().to_rgb
    clrs = [c(i) for i in clist]
    rvb = mcolors.LinearSegmentedColormap.from_list("", clrs)
    return rvb