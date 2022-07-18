# filter-regions

Filter genomic regions by score and proximity

## Installation

```bash
pip install .
```

or install with `Poetry`

```bash
poetry add filter-regions
```

## Usage

This module may be used as a library integrated with other Python projects, or as a standalone command-line tool. Here we describe its use as a Python library. Please run `filter-regions --help` for usage options on the command line.

### Input

There are two options for input. Input can be a file or Numpy vector (1d array) containing tab-delimited floating-point values on one line (input type of `vector`). Alternatively, input may be a BedGraph-formatted (BED3+1) file or Numpy matrix (2d array), where the fourth column contains floats (`bedgraph`).

### Output

Depending on the chosen input type (`vector` or `bedgraph`), by default, the output will be either a three- or four-column Pandas dataframe, respectively.

You can specify a ceiling on the number of elements reported via the `max_elements` parameter. If the number of elements found is greater than this property, the highest-scoring `max_elements` number of elements are reported. If there are fewer elements found, that smaller number is reported.

#### Vector

In the case of the `vector` input type, the output will contain half-open start and end indices (`[start, end)`) and the aggregated score over the half-open interval.

The difference of these indices will equal the `window_bins` parameter (e.g., 125 in current usage cases).

Here are the column names for the output dataframe:

 1. `Start`
 2. `Stop`
 3. `Score`

#### BedGraph

In the case of `bedgraph` input, the output will contain three columns for the chromosome, and start and stop positions (also using BED-like, half-open indexing). The fourth column will contain the aggregated score over the interval. 

The difference of the start and stop positions will equal the `window_bins` parameter, times the size of an individual bin (in units of nt). 

Here are the column names for the output dataframe:

 1. `Chromosome`
 2. `Start`
 3. `Stop`
 4. `Score`

Additional columns may be enabled via the `preserve_cols` parameter. Please read below for more detail.

### Filter method

There are three filter methods available: `pq` ([priority-queue](https://en.wikipedia.org/wiki/Priority_queue)), `wis` ([weighted-interval scheduling](https://en.wikipedia.org/wiki/Interval_scheduling#Weighted)), and `maxmean` (max-mean sweep).

The max-mean sweep method (MM) is a modification of the priority-queue method, using both the maximum and mean scores over each *window* to prioritize window selection. 

The primary difference between MM and PQ methods is that MM will use an aggregated score over the window to prioritize its selection, while PQ will only use the score at the centermost position within the window.

### Aggregation method

Once an interval has been selected to be filtered using whatever filter method, its assigned score value can be its original score — whatever the score is at the index in the vector, or bin in the bedgraph file — or an *aggregated* value that is calculated from applying a function on scores over the filtered window. 

The following functions are available via the `aggregation_method` parameter:

 * `min`
 * `max`
 * `mean`
 * `sum`
 * `median`
 * `variance`
 * `percentile`

The default is `max`, *i.e.*, given a filtered window, the score returned is the maximum value of the window. The other functions work accordingly, per their name. 

In the case of the `percentile` aggregation method, the `percentile` parameter may be specified as a value between 0 and 1. Its default is 0.95; however, another value may be specified. Using 0.5, for example, would be functionally equivalent to applying the `median` function.

### Additional columns

For either input type, if the `preserve_cols` flag is set to `True`, additional columns are reported in output:

 1. `OriginalIdx`
 2. `RollingMin`  
 3. `RollingMax` 
 4. `RollingMean`  
 5. `RollingSum`  
 6. `RollingMedian`  
 7. `RollingVariance`  
 8. `RollingPercentile`  
 9. `MethodIdx`

The `Rolling*` statistics are the results of applying the aggregation function to signal within each window. 

The `OriginalIdx` and `MethodIdx` values are indices indicating the starting bin index and the index used for window selection. This can be useful for evaluating which bin is used for selection.

#### Example

Here is a minimal example of how to use the installed package as a library in a standalone script, which reads in and filters a tab-delimited string of signal values from a regular file:

```
#!/usr/bin/env python

import filter_regions as fr

m = 'maxmean'
i = 'tests/test_example/scores.txt'
t = 'vector'
w = 125

f = fr.Filter(method=m,
              input=i,
              input_type=t,
              window_bins=w)
f.read()
f.filter()
o = f.output_df
print(o.head())
```

The `Filter` class property `output_df` specifies a Pandas dataframe on which all the usual methods may be called (*e.g.*, `head()` etc.). One can instead use the `write()` class method to send results to standard output or to a standard file, for example:

```
f.write() # send columns to standard output stream
```

Or:

```
f.write(output='/path/to/output') # write columns to a tab-delimited file
```

#### Example (advanced)

Other options may be used, *e.g.*, to pass in a BedGraph file, use the `mean` aggregation function, and report all window statistics:

```
#!/usr/bin/env python

import filter_regions as fr

m = 'maxmean'
i = 'tests/test_example/scores.bed'
t = 'bedgraph'
w = 125
a = 'mean'
p = True

f = fr.Filter(method=m, 
              input=i,
              input_type=t,
              window_bins=w,
              aggregation_method=a,
              preserve_cols=p)
f.read()
f.filter()
o = f.output_df
# ...
```

#### Example (Numpy, 1d)

Here is an example where a 1d array or vector of floats may be passed in:

```
#!/usr/bin/env python

import sys
import numpy as np
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
```

#### Example (Numpy, 2d)

Finally, here is an example of passing in a BedGraph-like 2d array or matrix:

```
#!/usr/bin/env python

import sys
import numpy as np
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
```
