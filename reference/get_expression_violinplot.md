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
  by = "group",
  value_transform = c("log2", "zscore", "none"),
  summarise = FALSE,
  facet_by = c("auto", "gene", "none"),
  sample_order = c("input", "group", "expression")
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

- by:

  Plot unit. Violin plots currently support only `"group"`, because a
  violin needs replicate-level distributions within groups rather than
  single values per sample.

- value_transform:

  One of `"log2"`, `"zscore"`, or `"none"`.

- summarise:

  Logical retained for compatibility. Violin plots always use
  replicate-level values, so `summarise = TRUE` is ignored with a
  warning.

- facet_by:

  Faceting mode: `"auto"` (default), `"gene"`, or `"none"`.

- sample_order:

  Ordering for sample-level x-axis display: `"input"`, `"group"`, or
  `"expression"` before values are grouped into violins.

## Value

A `ggplot2` object.

## Examples

``` r
v <- example_vista()
#> estimating size factors
#> estimating dispersions
#> gene-wise dispersion estimates
#> mean-dispersion relationship
#> final dispersion estimates
#> fitting model and testing
genes <- head(rownames(v), 4)
p <- get_expression_violinplot(v, genes = genes)
print(p)
```
