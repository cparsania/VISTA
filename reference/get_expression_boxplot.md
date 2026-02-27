# Plot gene expression distributions as boxplots

Displays per-sample or per-group distributions for selected genes using
normalized counts. When multiple genes are supplied and `facet = TRUE`,
facets are per gene (not per group). When `facet = FALSE` with multiple
genes, genes are shown on the x-axis and facets are per group.

## Usage

``` r
get_expression_boxplot(
  x,
  genes = NULL,
  sample_group = NULL,
  group_column = NULL,
  facet = TRUE,
  log_transform = TRUE,
  display_id = NULL,
  display_from = NULL,
  display_orgdb = NULL,
  facet_scales = "free_y",
  stats_group = FALSE,
  p.label = "p.signif",
  comparisons = NULL,
  pool_genes = FALSE,
  x_by = "group",
  facet_by = "none",
  fill_by = NULL
)
```

## Arguments

- x:

  A `VISTA` object.

- genes:

  Optional character vector of genes to display (≤20). Defaults to all
  genes.

- sample_group:

  Optional character vector specifying which groups (as defined by
  `group_column`) to include.

- group_column:

  Optional column name in `sample_info` used as the grouping variable.

- facet:

  Logical; facet the plot when `facet_by` is not `"none"`.

- log_transform:

  Logical; apply log2(x + 1) transform before plotting.

- display_id:

  Optional column in `rowData(x)` to use for gene labels (facets).

- display_from:

  Optional source ID type for mapping (used when `display_id` is not
  found in `rowData`).

- display_orgdb:

  Optional `OrgDb` object used for ID mapping when `display_id` is set
  but not found in `rowData`.

- facet_scales:

  Facet scales argument passed to `facet_wrap()` (default `"free_y"`).

- stats_group:

  Logical; add statistical comparisons between groups when `TRUE`.

- p.label:

  Label format for
  [`ggpubr::stat_compare_means()`](https://rpkgs.datanovia.com/ggpubr/reference/stat_compare_means.html).

- comparisons:

  Optional list of specific group comparisons for
  `stat_compare_means()`.

- pool_genes:

  Logical; when `TRUE`, pool selected genes into one distribution per
  x-axis category (scenario 1).

- x_by:

  When `pool_genes = TRUE`, either `"group"` or `"sample"` (x-axis and
  optional fill). When `pool_genes = FALSE`, either `"group"` or
  `"gene"` (x-axis for per-gene distributions).

- facet_by:

  Faceting control: for pooled genes, `"group"` or `"none"`; for
  per-gene, `"none"` (default) or `"gene"`.

- fill_by:

  When `pool_genes = TRUE`, either `"x"` (default) or `"group"` to force
  group colors even if x_by = "sample". When `pool_genes = FALSE`,
  either `"gene"` or `"group"`.

## Value

An object returned by this function.

## Examples

``` r
NULL
#> NULL
```
