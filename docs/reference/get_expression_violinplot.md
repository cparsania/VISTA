# Violin plot of expression values

Shows per-sample (or per-group) expression distributions as violins with
optional faceting by group.

## Usage

``` r
get_expression_violinplot(
  x,
  genes = NULL,
  sample_group = NULL,
  group_column = NULL,
  value_transform = c("log2", "zscore", "none"),
  summarise = FALSE,
  facet = TRUE
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

## Value

A `ggplot2` object.
