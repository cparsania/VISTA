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
  value_transform = c("log2", "zscore", "none"),
  summarise_groups = FALSE,
  facet_by_group = FALSE,
  color_palette = "Dark 3"
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

- value_transform:

  Transformation applied to expression values; one of `"log2"`,
  `"zscore"`, or `"none"`.

- summarise_groups:

  Logical; if `TRUE`, averages replicates per group before plotting.

- facet_by_group:

  Logical; facet the plot by the grouping column when `TRUE`.

- color_palette:

  Qualitative palette name used for gene colors.
