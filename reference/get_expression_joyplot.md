# Plot expression distributions as ridgelines

Shows per-group (or per-sample) expression distributions pooled across
the selected genes using ridge (joy) plots. Genes are pooled; no
faceting to keep shapes comparable.

## Usage

``` r
get_expression_joyplot(
  x,
  genes = NULL,
  sample_group = NULL,
  group_column = NULL,
  log_transform = TRUE,
  alpha = 0.7,
  scale = 1.2,
  y_by = c("group", "sample"),
  color_by = c("group", "sample")
)
```

## Arguments

- x:

  A `VISTA` object.

- genes:

  Optional character vector of genes to display. Defaults to all genes.

- sample_group:

  Optional character vector specifying which groups (as defined by
  `group_column`) to include.

- group_column:

  Optional column name in `sample_info` used as the grouping variable.

- log_transform:

  Logical; apply log2(x + 1) transform before plotting.

- alpha:

  Numeric transparency for fills.

- scale:

  Numeric scaling factor for ridges (passed to `geom_density_ridges()`).

- y_by:

  Either `"group"` (default) or `"sample"` to choose the y-axis (ridge)
  grouping.

- color_by:

  Either `"group"` (default) or `"sample"` to choose fill colors.

## Value

A `ggplot2` object.

## Examples

``` r
NULL
#> NULL
```
