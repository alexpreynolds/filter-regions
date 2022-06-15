# filter-regions

Filter genomic regions by score and proximity

## Installation

```bash
pip install -U filter-regions
```

or install with `Poetry`

```bash
poetry add filter-regions
```

Then you can run

```bash
filter-regions --help
```

or with `Poetry`:

```bash
poetry run filter-regions --help
```

## Usage

### Command-line

Example:

```
filter-regions --method maxmean --input input.bed > output.bed
```

The example above uses the max-mean sweep method to locate high-priority intervals. Replace `maxmean` with `pq` or `wis` to use the priority-queue or weighted-interval scheduling method for filtering.

For methods `pq` and `maxmean`, the `input.bed` file should be BED3+, where fourth and additional columns are a per-interval score vector. For the `wis` method, the `input.bed` file should be BedGraph or BED3+1; the fourth column should be a score value.

Depending on the particulars of your input, you can also override bin size and exclusion size:

```
filter-regions --method maxmean --input input.bed --bin-size 200 --exclusion-size 24800 > output.bed
```

In this example, the entire exclusion space is therefore 25kb (200 + 24800 nt). Filtered elements will not overlap within this space.

### Library 

Here is one example of how to use the installed package as a library in a standalone script:

```
#!/usr/bin/env python

import filter_regions as fr

m = 'maxmean'
i = '/path/to/scores.bed'
w = 125

f = fr.Filter(method=m, input=i, window_bins=w)
f.read()
f.filter()
o = f.output_df
print(o.head())
```

The `Filter` class property `output_df` returns a Pandas dataframe, on which all the usual methods may be called (e.g., `head()` etc.).

One can instead use the `write()` class method to send results to standard output or to a standard file:

```
f.write() # send to standard output stream
```

Or:

```
f.write(output='/path/to/output') # write to a file
```

#### Preserve columns

For the `pq` and `maxmean` methods, you may like to preserve additional columns in the output. Run the following after `read()` and before `filter()` to do so:

```
f.preserve_cols = True
```

Then `f.filter()` etc.