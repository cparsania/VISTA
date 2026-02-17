# Get Genes by Regulation

Extract gene IDs by regulation class from selected comparisons in a
VISTA object.

## Usage

``` r
get_genes_by_regulation(x, sample_comparisons, regulation = "Both")
```

## Arguments

- x:

  A VISTA object.

- sample_comparisons:

  Character vector of comparison names to include.

- regulation:

  One of "Up", "Down", "Both", or "All" (default: "Both").

## Value

A named list of character vectors (one per comparison).
