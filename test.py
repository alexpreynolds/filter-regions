#!/usr/bin/env python

import filter_regions as fr

m = 'pq'
# i = 'tests/test_example/scores.txt'
# t = 'vector'
i = 'tests/test_example/scores.10k.bed'
t = 'bedgraph'
w = 125
a = 'percentile'
c = 0.5
p = True

f = fr.Filter(method=m, input=i, input_type=t, aggregation_method=a, percentile=c, window_bins=w, preserve_cols=p)
# f = fr.Filter(method=m, input=i, input_type=t, window_bins=w)
f.read()
f.filter()
f.write(output=None)
