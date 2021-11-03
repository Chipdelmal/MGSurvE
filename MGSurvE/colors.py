#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap

def colorPaletteFromHexList(clist):
    """Generates a matplotlib-compliant cmap from a list of hex colors.

    Args:
        clist (list): List of hex codes (eg. '#f72585') to blend in the map.

    Returns:
        cmap: Matplotlib's colormap object.
    """
    c = mcolors.ColorConverter().to_rgb
    clrs = [c(i) for i in clist]
    rvb = mcolors.LinearSegmentedColormap.from_list("", clrs)
    return rvb