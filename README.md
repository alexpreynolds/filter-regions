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

It can be useful to try out the script with the example intervals:

```
unstarch tests/test_example/scores.bed.starch > tests/test_example/scores.bed
```

Then:

```
filter-regions --method maxmean --input tests/test_example/scores.bed > /tmp/scores.filtered_via_maxmean.bed
```

This uses the max-mean sweep method to locate high-priority intervals. Replace `maxmean` with `pq` or `wis` to use the priority-queue or weighted-interval scheduling method for filtering.

Depending on the particulars of your input, you can also override bin size and exclusion size:

```
filter-regions --method maxmean --input tests/test_example/scores.bed --bin-size 200 --exclusion-size 24800 > /tmp/scores.filtered_via_maxmean.bed
```

In this example, the entire exclusion space is therefore 25kb (200 + 24800 nt). Filtered elements will not overlap within this space.