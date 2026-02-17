# Donut chart of DEG counts (Up/Down) across comparisons

Donut chart of DEG counts (Up/Down) across comparisons

## Usage

``` r
get_deg_count_donutplot(
  x,
  sample_comparisons = NULL,
  show_counts = TRUE,
  show_percent = TRUE,
  percent_digits = 1,
  font_size = 12,
  fill_colors = c(Up = "red4", Down = "blue4"),
  facet_by = c("sample_comparisons", "none"),
  ncol = NULL
)
```

## Arguments

- x:

  A `VISTA` object containing DEG summaries.

- sample_comparisons:

  Optional character vector of comparison names to display.

- show_counts:

  Logical; include raw count labels in slices.

- show_percent:

  Logical; include percent labels in slices.

- percent_digits:

  Integer number of decimals used for percentage labels.

- font_size:

  Numeric base font size for the plot.

- fill_colors:

  Named vector giving fill colors for `"Up"` and `"Down"` slices.

- facet_by:

  Either `"sample_comparisons"` (default) to draw one donut per
  comparison, or `"none"` for a single donut.

- ncol:

  Optional number of columns when faceting.
