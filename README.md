## Filter HTSeq Counts Matrices

The `filter-counts` tool filters the genes in an HTSeq counts matrix file.

## Options & Arguments

filter-counts is called as:

~~~
USAGE:
    filter-counts [FLAGS] [OPTIONS] <path>

FLAGS:
    -i, --filter-identical    Filter out genes with zero variance (i.e. with all values identical)
    -h, --help                Prints help information
    -V, --version             Prints version information
    -v, --verbose             Provide verbose output. supply multiple times to increase verbosity

OPTIONS:
    -z, --max-zerocount <n>    Maximum number of zero counts permitted in a single gene
    -m, --min-count <n>        Minimum total gene count

ARGS:
    <path>    Input counts file
~~~

## Filters

Each gene if filtered based on the following criteria in order:

1. The total counts for the gene must be ≥ the value specified by `-m` (if specified);
2. The total number of zero counts for the gene must be ≤ the value specified by `-z`;
3. The gene variance must be > 0 (if `-i` is specified).

## Licence
These tools are released under the [MIT License](https://opensource.org/licenses/MIT).

## Building

Before installation, you'll need to install [Rust](https://www.rust-lang.org/).

~~~bash
$ git clone https://github.com/alastair-droop/filter-counts.git
$ cd filter-counts
$ cargo install --path .
~~~
