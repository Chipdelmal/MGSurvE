#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

distances = [4.572, 9.144, 18.288, 36.576, 54.864, 73.152]
traps = {
    'An. melas': {
        'calves':   [9166, 2712, 1419, 591, 181,  97],
        'co2':      [1245,  505,  195,  83,  54,  26],
        'control':  [ 200,  182,  217, 185, 164, 131]
    },
    'An. other': {
        'calves':   [ 618,  110,   40,   8,   6,   0],
        'co2':      [ 103,   25,   11,   4,   0,   1],
        'control':  [   4,    5,    3,   3,   1,   2]
    },
    'Ae. other': {
        'calves':   [1139,  292,  92,   31,  19,  12],
        'co2':      [ 238,   75,  39,    7,  12,  12],
        'control':  [  10,    9,  16,   14,  10,  18]
    }
}

(fig, ax) = plt.subplots(figsize=(10, 5))
for typ in ['calves', 'co2', 'control']:
    plt.plot(
        distances, 
        np.array(traps['An. other'][typ])+np.array(traps['An. melas'][typ])
    )
ax.set_xlim(0, distances[-1])
ax.set_ylim(0, 10250)


(fig, ax) = plt.subplots(figsize=(10, 5))
for typ in ['calves', 'co2', 'control']:
    plt.plot(
        distances, 
        np.array(traps['Ae. other'][typ])
    )
ax.set_xlim(0, distances[-1])
ax.set_ylim(0, 1500)