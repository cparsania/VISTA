# Gene expression line plot

Plots normalized expression for selected genes across samples or
summarized groups with optional transformations and group faceting.

## Usage

``` r
get_expression_lineplot(
  x,
  genes,
  sample_group = NULL,
  group_column = NULL,
  by = c("sample", "group"),
  value_transform = c("log2", "zscore", "none"),
  facet_by = c("none", "group", "gene"),
  sample_order = c("input", "group", "expression"),
  palette = NULL,
  colors = NULL,
  line_width = 1,
  point_size = 2,
  base_size = 12
)
```

## Arguments

- x:

  A `VISTA` object.

- genes:

  Character vector of gene identifiers to plot.

- sample_group:

  Optional character vector specifying which groups (values taken from
  `group_column`) to include.

- group_column:

  Optional column name in `sample_info` defining the grouping/faceting
  variable.

- by:

  Plot unit: `"sample"` (default) or `"group"` to average replicates
  before plotting.

- value_transform:

  Transformation applied to expression values; one of `"log2"`,
  `"zscore"`, or `"none"`.

- facet_by:

  Faceting mode: `"none"` (default), `"group"`, or `"gene"`.

- sample_order:

  Ordering used for sample-level plots: `"input"`, `"group"`, or
  `"expression"`.

- palette:

  Optional qualitative palette name used for gene colours.

- colors:

  Optional named character vector of manual gene colours.

- line_width:

  Line width.

- point_size:

  Point size.

- base_size:

  Base theme size.

## Value

An object returned by this function.

## Examples

``` r
v <- example_vista()
genes <- head(rownames(v), 3)
p <- get_expression_lineplot(v, genes = genes)
print(p)
```
