# Plot pairwise correlations between samples

Uses GGally::ggpairs on normalized expression to display correlations
among samples from selected groups/genes.

## Usage

``` r
get_pairwise_corr_plot(
  x,
  sample_group = NULL,
  group_column = NULL,
  genes = NULL,
  sample_order = c("input", "group"),
  value_transform = c("log2", "none"),
  title = "Pairwise Sample Correlations"
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

- sample_order:

  Ordering for selected samples: `"input"` keeps the current order,
  while `"group"` sorts by `group_column`.

- value_transform:

  Value transformation: `"log2"` (default) or `"none"`.

- title:

  Plot title.

## Value

An object returned by this function.

## Examples

``` r
NULL
#> NULL
```
