# type: ignore[attr-defined]

import bisect
import collections
import errno
import os
import sys
import timeit
from enum import Enum, EnumMeta, unique
from io import StringIO

import natsort as ns
import numpy as np
import pandas as pd
import pyranges as pr
import typer
from rich.console import Console


APP_NAME = 'filter-regions'
MAX_ELEMENTS = 100_000


class FilterMethodMeta(EnumMeta): 
    def __contains__(cls, item): 
        return item in [v.value for v in cls.__members__.values()]


@unique
class FilterMethod(str, Enum, metaclass=FilterMethodMeta):
    pq = "pq"
    wis = "wis"
    maxmean = "maxmean"


class FilterMethodDescription(str, Enum, metaclass=FilterMethodMeta):
    pq = "Priority-Queue (PQ)"
    wis = "Weighted Interval Scheduler (WIS)"
    maxmean = "Max-Mean Sweep (MaxMean)"


class Filter:
    def __init__(self, method, input, window_bins=None, bin_size=200, exclusion_size=24800, preserve_cols=False) -> None:
        self.method: str = method
        self.input: str = input
        self.window_bins: int = window_bins if window_bins else None
        self.bin_size: int = bin_size
        self.exclusion_size: int = exclusion_size
        self.input_df: pd.DataFrame = None
        self.output_df: pd.DataFrame = None
        self.preserve_cols: bool = preserve_cols

    def input_df(self, df: pd.DataFrame) -> None:
        self.input_df = df

    def output_df(self, df: pd.DataFrame) -> None:
        self.output_df = df

    def preserve_cols(self, flag: bool) -> None:
        self.preserve_cols = flag

    def read(self) -> None:
        df = None
        console.print(f"[bold blue]{APP_NAME}[/] [bold green]Run[/] Reading input file into dataframe...")
        if self.method == 'pq' or self.method == 'maxmean':
            df = pd.read_csv(self.input, sep='\t', header=None, names=['Chromosome', 'Start', 'End', 'Score'])
        elif self.method == 'wis':
            df = pr.read_bed(self.input)
        self.input_df = df
        self.window_bins = self.window_bins if self.window_bins else (self.bin_size + self.exclusion_size) // self.bin_size # window size (in units of bins, not nt)

    def write(self, output: None) -> None:
        o = StringIO()
        try:
            self.output_df.to_csv(o, sep="\t", index=False, header=False)
            if not output:
                sys.stdout.write(f"{o.getvalue()}")
            else:
                with open(output, 'wb') as ofh:
                    ofh.write(f"{o.getvalue()}".encode("utf-8"))
        except AttributeError as err:
            console.print(f"[bold blue]{APP_NAME}[/] [bold red]Error[/] Could not write output \[{err}]")
            sys.exit(errno.EINVAL)
    
    def output_df(self) -> pd.DataFrame:
        return self.output_df

    def filter(self) -> None:
        try:
            console.print(f"[bold blue]{APP_NAME}[/] [bold green]Run[/] Filtering with {self.method} method...")
            start = timeit.default_timer()
            run = getattr(self, self.method)
            self.output_df = run()
            end = timeit.default_timer()
            console.print(f"[bold blue]{APP_NAME}[/] [bold green]Run[/] Method completed in: {(end - start):.2f} s")
            console.print(f"[bold blue]{APP_NAME}[/] [bold green]Run[/] Method found: {len(self.output_df.index)} elements")
        except TypeError:
            console.print(f"[bold blue]{APP_NAME}[/] [bold red]Error[/] Method could not be applied")
            sys.exit(errno.EINVAL)
    
    def pq(self) -> pd.DataFrame:
        maxOnly = True
        return self.maxmean(maxOnly)

    def wis(self) -> pd.DataFrame:
        d = self.input_df
        m = d.merge()
        df = d.as_df()
        df = df.rename(columns={"Name" : "Score"})
        # update the start and end loci based on window size (bin units)
        df["Start"] = df.Start.shift(self.window_bins // 2)
        df["End"] = df.End.shift(-(self.window_bins // 2)) if self.window_bins % 2 else df.End.shift(-(self.window_bins // 2 - 1))
        df = df.dropna()
        df = df[df['Start'] < df['End']]
        df["Start"] = df["Start"].astype(int)
        df["End"] = df["End"].astype(int)
        # intervals in the df must be sorted for the trace to work correctly (natural sort on chromosome)
        df = df.sort_values(by=["Chromosome", "Start", "End"])
        df = df.reindex(index=ns.order_by_index(df.index, ns.index_natsorted(zip(df.Chromosome, df.Start, df.End))))
        df = df.reset_index(drop=True)
        n = len(df.index)
        # build accumulator-based dict for calculation of absolute coordinates
        acc_chr = m.as_df().loc[:, "Chromosome"].copy().values
        acc_abs = np.concatenate([[0], m.as_df().loc[:, 'End'].copy().values[:-1]])
        j = len(acc_abs) - 1
        while j > 0:
            acc_abs[j] = np.sum(acc_abs[0:j+1])
            j -= 1
        acc = dict(zip(acc_chr, acc_abs))
        # for every interval j, compute the rightmost disjoint interval i, where i < j.
        s = df.loc[:, "Start"].copy().values
        e = df.loc[:, "End"].copy().values
        # translate per-chr coord space to abs coord space
        c = df.loc[:, "Chromosome"].copy().values
        translate_genomic_coords_to_abs = np.vectorize(lambda x, y: y + acc[x])
        s = translate_genomic_coords_to_abs(c, s)
        e = translate_genomic_coords_to_abs(c, e)
        # look for nearest disjoint interval to left of current interval
        p = []
        for j in range(n):
            i = bisect.bisect_right(e, s[j]) - 1
            p.append(i)
        # set up initial opt(imum) table values.
        opt = collections.defaultdict(int)
        opt[-1] = 0
        opt[0] = 0
        # forwards pass to get best-scoring path
        for j in range(1, n):
            score = df.loc[df.index[j], "Score"]
            opt[j] = max(score + opt[p[j]], opt[j - 1])
        # backwards trace to retrieve path
        q = []
        j = n - 1
        print('{}'.format(j))
        while j >= 0:
            score = df.loc[df.index[j], "Score"]
            if score + opt[p[j]] > opt[j - 1]:
                q.append(j)
                j = p[j] # jump to the nearest disjoint interval to the left
            else:
                j -= 1 # try the "next" interval, one to the left
        # sort indices of qualifiers
        q.sort()
        r = df.iloc[q]
        return r
    
    def maxmean(self, maxOnly: bool=False) -> pd.DataFrame:
        df = self.input_df.iloc[:, :3]
        df[3] = self.input_df.iloc[:, 3:].sum(axis=1)
        df.columns = ["Chromosome", "Start", "End", "Score"]
        # save original indices
        df["ScoresIdx"] = df.index
        # update the start and end loci based on window size (bin units)
        df["Start"] = df.Start.shift(self.window_bins // 2)
        df["End"] = df.End.shift(-(self.window_bins // 2)) if self.window_bins % 2 else df.End.shift(-(self.window_bins // 2 - 1))
        df = df.dropna()
        df["Start"] = df["Start"].astype(int)
        df["End"] = df["End"].astype(int)
        df = df.reset_index(drop=True)
        # get rolling means and maxes
        df['Max'] = df['Score'].rolling(self.window_bins, center=True).max()
        df['Mean'] = df['Score'].rolling(self.window_bins, center=True).mean()
        # get rid of edges
        df = df.dropna().reset_index(drop=True)
        # drop regions which fall over two chromosomes
        df = df.drop(np.where(df.iloc[:, 1].to_numpy() >= df.iloc[:, 2].to_numpy())[0]).reset_index(drop=True)
        # save new indices before sorting
        df['MethodIdx'] = df.index
        # sort by max, by mean, and then by score
        df = df.sort_values(by=["Max", "Mean", "Score"], ascending=False) if not maxOnly else df.sort_values(by=["Max", "Score"], ascending=False)
        # run max-mean sweep algorithm
        n = len(df)
        hits = np.zeros(n, dtype=bool)
        indices = []
        k = MAX_ELEMENTS
        for loc in range(n):
            if k <= 0:
                break
            j = df.iloc[loc, 7]
            start = (j - self.window_bins // 2) if (j - self.window_bins // 2) > 0 else 0
            if self.window_bins % 2: # if there are an odd number of bins, have to add one to stop                                                                                              
                stop = (j + self.window_bins // 2 + 1) if (j + self.window_bins // 2 + 1) <= n else n
            else:
                stop = (j + self.window_bins // 2) if (j + self.window_bins // 2) <= n else n
            if not np.any(hits[start:stop]):
                hits[start:stop] = True
                indices.append(j)
                k -= 1
        df = df.loc[indices]
        if not maxOnly: df.Score = df.Max
        df = df.sort_values(by=["ScoresIdx"], ascending=True)
        if not self.preserve_cols: df = df[["Chromosome", "Start", "End", "Score"]]
        df = df.loc[df["Score"] > 0]
        df = df.reset_index(drop=True)
        return df

app = typer.Typer(
    name="filter-regions",
    help="Filter genomic regions by score and proximity",
    add_completion=False,
)

# write console messages to standard error
console = Console(stderr=True)


@app.command(name="")
def main(
    method: str = typer.Option(
        ...,
        "-m",
        "--method",
        case_sensitive=True,
        help="Filter method (pq|wis|maxmean)"
    ),
    input: str = typer.Option(
        ...,
        "-i",
        "--input",
        case_sensitive=True,
        help="Input filename path"
    ),
    window_bins: int = typer.Option(
        None,
        "-w",
        "--window-bins",
        help="Window of bins that excludes overlap (bins)"
    ),
    bin_size: int = typer.Option(
        200,
        "-b",
        "--bin-size",
        help="Bin size (nt)"
    ),
    exclusion_size: int = typer.Option(
        24800,
        "-x",
        "--exclusion-size",
        help="Exclusion size, up- and downstream of bin (nt)"
    ),
    preserve_cols: bool = typer.Option(
        False,
        "-p",
        "--preserve-cols",
        help="Keep all columns"
    ),
) -> None:
    """Filter regions by specified method."""
    
    if method not in FilterMethod:
        console.print(f"[bold blue]{APP_NAME}[/] [bold red]Error[/] Must specify valid selection method \[{'|'.join([e.value for e in FilterMethod])}]")
        sys.exit(errno.EINVAL)

    if not os.path.exists(input):
        console.print(f"[bold blue]{APP_NAME}[/] [bold red]Error[/] Must specify valid input path \[{input}]")
        sys.exit(errno.EINVAL)

    f = Filter(method, input, window_bins, bin_size, exclusion_size, preserve_cols)
    f.read()
    f.filter()
    f.write(output=None)

if __name__ == "__main__":
    app()
