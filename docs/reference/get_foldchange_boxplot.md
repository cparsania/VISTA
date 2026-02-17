# Plot fold-change distributions across comparisons

Builds boxplots of log2 fold changes for selected genes and comparisons,
optionally adding statistics.

## Usage

``` r
get_foldchange_boxplot(
  x,
  genes = NULL,
  sample_comparisons = NULL,
  sample_comparison = NULL,
  facet = TRUE,
  p.label = "p.signif",
  stats_group = FALSE
)
```

## Arguments

- x:

  A `VISTA` object containing differential expression results.

- genes:

  Optional character vector of gene IDs to include.

- sample_comparisons:

  Optional character vector of comparison names to plot.

- sample_comparison:

  Deprecated alias for `sample_comparisons`.

- facet:

  Logical; facet each comparison when `TRUE`.

- p.label:

  Label type passed to
  [`ggpubr::stat_compare_means()`](https://rpkgs.datanovia.com/ggpubr/reference/stat_compare_means.html).

- stats_group:

  Logical; add pairwise statistical tests when `TRUE`.
