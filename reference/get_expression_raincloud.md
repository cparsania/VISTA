# Raincloud plot of expression values

Combines a violin, boxplot, and jittered points per sample/group to show
distribution, summary, and individual values.

## Usage

``` r
get_expression_raincloud(
  x,
  genes = NULL,
  sample_group = NULL,
  group_column = NULL,
  value_transform = c("log2", "zscore", "none"),
  summarise = FALSE,
  facet = TRUE,
  point_alpha = 0.5,
  point_size = 1.5
)
```

## Arguments

- x:

  A `VISTA` object.

- genes:

  Optional character vector of gene IDs to include; defaults to all.

- sample_group:

  Optional subset of groups (values of `group_column`) to keep.

- group_column:

  Grouping column in `sample_info`; defaults to the stored grouping.

- value_transform:

  One of `"log2"`, `"zscore"`, or `"none"`.

- summarise:

  Logical; if `TRUE`, averages replicates per group before plotting.

- facet:

  Logical; facet by group.

- point_alpha:

  Alpha for jittered points.

- point_size:

  Point size for jittered points.

## Value

A `ggplot2` object.
