#!/usr/bin/env python

import sys

import numpy as np

sys.path.insert(1, '../')
import filter_regions as fr

m = 'wis'
i = np.array([
    ['chr1', 0, 100, 1.0],
    ['chr1', 100, 200, 2.0],
    ['chr1', 200, 300, 5.0],
    ['chr1', 300, 400, 10.0],
    ['chr1', 400, 500, 5.0],
    ['chr1', 500, 600, 5.0],
    ['chr1', 600, 700, 4.0],
    ['chr1', 700, 800, 8.0],
    ['chr2', 0, 100, 2.0],
    ['chr2', 100, 200, 1.0],
    ['chr2', 200, 300, 10.0],
    ['chr2', 300, 400, 1.0],
    ['chr2', 400, 500, 1.0],
    ['chr2', 500, 600, 1.0],
    ['chr2', 600, 700, 5.0],
    ['chr2', 700, 800, 5.0]
])
t = 'bedgraph'
w = 3
a = 'max'
p = True
q = False

f = fr.Filter(method=m, input=i, input_type=t, aggregation_method=a, window_bins=w, preserve_cols=p, quiet=q)
f.read()
f.filter()
f.write(output=None)
