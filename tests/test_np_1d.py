#!/usr/bin/env python

import sys

import numpy as np

sys.path.insert(1, '../')
import filter_regions as fr

m = 'maxmean'
i = np.array([1.0, 2.0, 5.0, 10.0, 5.0, 5.0, 4.0, 8.0, 2.0, 1.0, 10.0, 1.0, 1.0, 1.0, 5.0])
t = 'vector'
w = 3
a = 'max'
p = True
q = False

f = fr.Filter(method=m,
              input=i,
              input_type=t,
              aggregation_method=a,
              window_bins=w,
              preserve_cols=p,
              quiet=q)
f.read()
f.filter()
f.write(output=None)
