# Plot fold-change distributions across comparisons

Builds boxplots of log2 fold changes for selected genes and comparisons,
optionally adding statistics.

## Usage

``` r
get_foldchange_boxplot(
  x,
  genes = NULL,
  sample_comparisons = NULL,
  facet_by = c("auto", "comparison", "none"),
  p.label = "p.signif",
  stats_group = FALSE,
  stats_method = "t.test"
)
```

## Arguments

- x:

  A `VISTA` object containing differential expression results.

- genes:

  Optional character vector of gene IDs to include.

- sample_comparisons:

  Optional character vector of comparison names to plot.

- facet_by:

  Faceting mode: `"auto"` (default), `"comparison"`, or `"none"`.

- p.label:

  Label type passed to
  [`ggpubr::stat_compare_means()`](https://rpkgs.datanovia.com/ggpubr/reference/stat_compare_means.html).

- stats_group:

  Logical; add pairwise statistical tests when `TRUE`.

- stats_method:

  Statistical method passed to
  [`ggpubr::stat_compare_means()`](https://rpkgs.datanovia.com/ggpubr/reference/stat_compare_means.html).

## Value

An object returned by this function.

## Examples

``` r
NULL
#> NULL
```
