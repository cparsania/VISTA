# Plot cell type proportions

Plot cell type proportions

## Usage

``` r
plot_celltype_barplot(x, group_column = NULL, samples = NULL, font_size = 12)
```

## Arguments

- x:

  A VISTA object.

- group_column:

  A column in colData(x) to group by (e.g., "condition").

- samples:

  Optional character vector of sample names to include.

- font_size:

  Base font size.

## Value

A ggplot object.
