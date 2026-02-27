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
  facet_clusters = FALSE,
  line_transparency = 0.5,
  show_average_line = TRUE,
  average_line_color = NULL,
  average_line_size = 1,
  average_line_summary_method = "median"
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

- facet_clusters:

  Logical; facet the plot by cluster when k-means clustering is
  requested.

- line_transparency:

  Numeric alpha applied to individual gene lines.

- show_average_line:

  Logical; overlay a summary line per cluster when `TRUE`.

- average_line_color:

  Color used for the summary line. When `NULL`, uses the first
  comparison color (if stored) for consistency across plots.

- average_line_size:

  Numeric line width for the summary line.

- average_line_summary_method:

  Character string selecting `"median"` or `"mean"` for the summary
  statistic.

## Value

An object returned by this function.

## Examples

``` r
NULL
#> NULL
```
