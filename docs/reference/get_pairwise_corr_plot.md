# Plot pairwise correlations between samples

Uses GGally::ggpairs on normalized expression to display correlations
among samples from selected groups/genes.

## Usage

``` r
get_pairwise_corr_plot(
  x,
  sample_group = NULL,
  group_column = NULL,
  genes = NULL
)
```

## Arguments

- x:

  A `VISTA` object.

- sample_group:

  Optional character vector of groups (from the column specified by
  `group_column`) used to subset samples.

- group_column:

  Optional column name in `sample_info` defining the grouping used for
  filtering.

- genes:

  Optional character vector of gene IDs to include; defaults to all
  genes.
