# Plot enrichment results using -log10(FDR)

Generates a dot plot of enrichment results for `enrichResult`,
`gseaResult`, or `compareClusterResult` objects (including those
returned by [`enrichMsigDB()`](enrichMsigDB.md)). Points are sized by
gene/set count and coloured by -log10(FDR). For compareCluster results,
the plot is faceted by cluster with top terms selected per cluster.

## Usage

``` r
get_enrichment_plot(x, top_n = 10, title = NULL)
```

## Arguments

- x:

  An object of class `enrichResult`, `gseaResult`, or
  `compareClusterResult`.

- top_n:

  Integer; number of top terms to plot (per cluster for compareCluster).

- title:

  Optional plot title.

## Value

A `ggplot2` object.
