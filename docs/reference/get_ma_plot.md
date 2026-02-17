# Generate MA plot from a VISTA object

Create an MA plot (log2 fold change vs mean expression) for a selected
comparison contained in a VISTA object. Genes are coloured by their
regulation class and the top results can be optionally labeled with gene
IDs.

## Usage

``` r
get_ma_plot(
  x,
  sample_comparison,
  point_size = 1.2,
  alpha = 0.6,
  fill_colors = c(Up = "#a40000", Down = "#16317d", Other = "gray70"),
  topn = 0,
  label_size = 3,
  repair = FALSE,
  display_id = NULL,
  display_from = NULL,
  display_orgdb = NULL
)
```

## Arguments

- x:

  A [VISTA](VISTA-class.md) object.

- sample_comparison:

  Character scalar naming the comparison to plot. Must match one of
  `names(comparisons(x))`.

- point_size:

  Numeric point size. Default: 1.2.

- alpha:

  Numeric transparency (0-1). Default: 0.6.

- fill_colors:

  Named character vector of colors for `"Up"`, `"Down"`, and `"Other"`
  genes.

- topn:

  Integer number of genes to label. Default: 0.

- label_size:

  Text size for labels. Default: 3.

- repair:

  Logical; if `TRUE`, attempt to shorten gene identifiers to symbols by
  stripping prefixes. Default: `FALSE`.

- display_id:

  Optional ID/column name to use for labels. If supplied and present in
  `rowData(x)`, those values are used; otherwise falls back to ID
  mapping.

- display_from:

  Optional source ID type for mapping (used when `display_id` is not
  found in `rowData`).

- display_orgdb:

  Optional `OrgDb` used for ID mapping when `display_id` is set but not
  found in `rowData`.

## Value

A [`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html) MA
plot.
