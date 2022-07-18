#!/usr/bin/env python

import sys

import numpy as np

sys.path.insert(1, '../')
import filter_regions as fr

m = 'maxmean'
i = np.array([
    ['chr1', 0, 100, 1.0],
    ['chr1', 100, 200, 2.0],
    ['chr1', 200, 300, 5.0],
    ['chr1', 300, 400, 10.0],
    ['chr1', 400, 500, 5.0],
    ['chr1', 500, 600, 5.0],
    ['chr1', 600, 700, 4.0],
    ['chr1', 700, 800, 8.0],
    ['chr1', 800, 900, 2.0],
    ['chr1', 900, 1000, 1.0],
    ['chr1', 1000, 1100, 10.0],
    ['chr1', 1100, 1200, 1.0],
    ['chr1', 1200, 1300, 1.0],
    ['chr1', 1300, 1400, 1.0],
    ['chr1', 1400, 1500, 5.0]
])
t = 'bedgraph'
w = 3
a = 'max'
p = True
q = False
e = 123

f = fr.Filter(method=m, 
              input=i, 
              input_type=t, 
              aggregation_method=a, 
              window_bins=w, 
              preserve_cols=p, 
              max_elements=e,
              quiet=q)
f.read()
f.filter()
f.write(output=None)
