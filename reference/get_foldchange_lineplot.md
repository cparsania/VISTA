# Fold-change line plot across comparisons

Plots log2 fold-change trajectories for selected genes across multiple
comparisons, optionally clustering genes.

## Usage

``` r
get_foldchange_lineplot(
  x,
  sample_comparisons,
  genes = NULL,
  km = NULL,
  facet_by = c("none", "cluster"),
  alpha = 0.5,
  show_summary = TRUE,
  summary_color = NULL,
  summary_linewidth = 1,
  summary_fun = c("median", "mean")
)
```

## Arguments

- x:

  A `VISTA` object containing differential expression results.

- sample_comparisons:

  Character vector of comparison names to include.

- genes:

  Optional character vector of gene identifiers to plot. Defaults to all
  genes.

- km:

  Optional integer specifying the number of k-means clusters to compute;
  `NULL` disables clustering.

- facet_by:

  Faceting mode: `"none"` (default) or `"cluster"` when k-means
  clustering is requested.

- alpha:

  Numeric alpha applied to individual gene lines.

- show_summary:

  Logical; overlay a summary line per cluster when `TRUE`.

- summary_color:

  Color used for the summary line. When `NULL`, uses the first
  comparison color (if stored) for consistency across plots.

- summary_linewidth:

  Numeric line width for the summary line.

- summary_fun:

  Character string selecting `"median"` or `"mean"` for the summary
  statistic.

## Value

An object returned by this function.

## Examples

``` r
NULL
#> NULL
```
