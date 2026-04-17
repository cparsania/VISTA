# Fold-change line plot across comparisons

Plots log2 fold-change trajectories for selected genes across multiple
comparisons, optionally clustering genes.

## Usage

``` r
get_foldchange_lineplot(
  x,
  sample_comparisons,
  genes = NULL,
  display_id = NULL,
  display_from = NULL,
  display_orgdb = NULL,
  km = NULL,
  facet_by = c("none", "cluster"),
  facet_scales = "fixed",
  facet_nrow = NULL,
  facet_ncol = NULL,
  alpha = 0.5,
  palette = NULL,
  show_summary = TRUE,
  summary_color = NULL,
  summary_linewidth = 1,
  summary_fun = c("median", "mean"),
  base_size = 14
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

- display_id:

  Optional ID/column name used to interpret `genes` and map the returned
  cluster table to display-friendly labels.

- display_from:

  Optional source ID type for mapping when `display_id` is not present
  in `rowData(x)`.

- display_orgdb:

  Optional `OrgDb` object used for identifier mapping when `display_id`
  is not present in `rowData(x)`.

- km:

  Optional integer specifying the number of k-means clusters to compute;
  `NULL` disables clustering.

- facet_by:

  Faceting mode: `"none"` (default) or `"cluster"` when k-means
  clustering is requested.

- facet_scales:

  Facet scales argument passed to `facet_wrap()` when
  `facet_by = "cluster"` (default `"fixed"`).

- facet_nrow, facet_ncol:

  Optional layout passed to `facet_wrap()` when faceting.

- alpha:

  Numeric alpha applied to individual gene lines.

- palette:

  Optional named or unnamed color vector used for cluster lines.

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

- base_size:

  Numeric base theme size.

## Value

An object returned by this function.

A list with `plot` (the `ggplot2` object) and `clustered_data`
(gene-to-cluster assignments).

## Examples

``` r
v <- example_vista()
comp <- names(comparisons(v))[1]
genes <- head(as.character(comparisons(v)[[comp]]$gene_id), 5)
p <- get_foldchange_lineplot(v, sample_comparison = comp, genes = genes)
print(p)
#> $plot
#> `geom_line()`: Each group consists of only one observation.
#> ℹ Do you need to adjust the group aesthetic?
#> `geom_line()`: Each group consists of only one observation.
#> ℹ Do you need to adjust the group aesthetic?

#> 
#> $clustered_data
#>                         gene_id    display_gene cluster
#> ENSG00000000003 ENSG00000000003 ENSG00000000003     All
#> ENSG00000000419 ENSG00000000419 ENSG00000000419     All
#> ENSG00000000457 ENSG00000000457 ENSG00000000457     All
#> ENSG00000000460 ENSG00000000460 ENSG00000000460     All
#> ENSG00000000971 ENSG00000000971 ENSG00000000971     All
#> 
```
