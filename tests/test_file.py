#!/usr/bin/env python

import sys

sys.path.insert(1, '../')
import filter_regions as fr

m = 'wis'
# i = 'test_example/scores.10k.txt'
# t = 'vector'
i = 'test_example/scores.bed'
t = 'bedgraph'
a = 'max'
w = 125
p = True
e = 125_000

f = fr.Filter(method=m,
              input=i,
              input_type=t,
              aggregation_method=a,
              window_bins=w,
              preserve_cols=p,
              max_elements=e)
f.read()
f.filter()
f.write(output=None)
