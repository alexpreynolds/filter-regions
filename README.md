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

This module may be used as a library integrated with other Python projects, or as a standalone binary. Here we describe its use as a Python library.

### Input

Input can be a file containing tab-delimited floating-point values on one line (input type of `vector`) or a BedGraph-formatted (BED3+) file, where the fourth column contains floats (`bedgraph`).

### Output

Depending on input type (`vector` or `bedgraph`), by default, the output will be either three- or four-column values, respectively. 

When used as a library, the output will be a Pandas dataframe with column names described below.

In the case of the `vector` input type, the output will contain a half-open start and end index and the aggregated score over that index (`[start, end)`). The difference of these indices will equal the `window_bins` parameter (e.g., 125 in current usage cases). Here are column names for this output:

 1. Start
 2. Stop
 3. Score

In the case of `bedgraph` input, the output will contain three columns for the chromosome, and start and stop positions (also half-open coordinate-indexed). The fourth column will contain the aggregated score over the interval. The difference of the start and stop positions will equal the `window_bins` parameter, times the size of an individual bin (units of nt). Here are column names for this type of output:

 1. Chromosome
 2. Start
 3. Stop
 4. Score

Additional columns can be enabled via the `preserve_cols` parameter; see below for more detail.

### Filter method

There are three filter methods available: `pq` ([priority-queue](https://en.wikipedia.org/wiki/Priority_queue)), `wis` ([weighted-interval scheduling](https://en.wikipedia.org/wiki/Interval_scheduling#Weighted)), and `maxmean` (max-mean sweep).

The max-mean sweep method (MM) is a modification of the priority-queue method, using both the maximum and mean scores over each *window* to prioritize window selection. 

The primary difference between MM and PQ methods is that MM will use this aggregated "score" over the window to prioritize its selection, while PQ will only use the score at the centermost position within the window.

### Aggregation method

Once an interval is selected to be filtered, its score value can be its original score (whatever the score is at the index in the vector, or bin in the bedgraph file), or an *aggregated* value calculated from applying a function on scores over the filtered window. The following functions are available via the `aggregation_method` parameter:

 * `min`
 * `max`
 * `mean`
 * `sum`
 * `median`
 * `variance`
 * `percentile`

The default is `max`, *i.e.*, given a filtered window, the score returned is the maximum value of the window. The other functions work accordingly, per their name. 

In the case of `percentile`, the `percentile` parameter may be specified as a value between 0 and 1. Its default is 0.95. However, another value may be specified; using 0.5, for example, would be equivalent to the `median` function.

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

#### Example

Here is a minimal example of how to use the installed package as a library in a standalone script:

```
#!/usr/bin/env python

import filter_regions as fr

m = 'maxmean'
i = 'tests/test_example/scores.txt'
t = 'vector'
w = 125

f = fr.Filter(method=m, input=i, input_type=t, window_bins=w)
f.read()
f.filter()
o = f.output_df
print(o.head())
```

The `Filter` class property `output_df` returns a Pandas dataframe, on which all the usual methods may be called (e.g., `head()` etc.). One can instead use the `write()` class method to send results to standard output or to a standard file, for example:

```
f.write() # send columns to standard output stream
```

Or:

```
f.write(output='/path/to/output') # write columns to a tab-delimited file
```

#### Example (advanced)

Other options can be used, e.g., to pass in a BedGraph file and use the `mean` aggregation function:

```
#!/usr/bin/env python

import filter_regions as fr

m = 'maxmean'
i = 'tests/test_example/scores.bed'
t = 'bedgraph'
w = 125
a = 'mean'
p = True

f = fr.Filter(method=m, input=i, input_type=t, window_bins=w, aggregation_method=a, preserve_cols=p)
f.read()
f.filter()
# ...
```