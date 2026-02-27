# Barplot of DEG counts (Up/Down) across comparisons

Barplot of DEG counts (Up/Down) across comparisons

## Usage

``` r
get_deg_count_barplot(
  x,
  sample_comparisons = NULL,
  show_counts = TRUE,
  font_size = 12,
  fill_colors = c(Up = "red4", Down = "blue4"),
  facet_by = c("none", "regulation", "sample_comparisons", "sample_comparison"),
  facet_scale = "fixed"
)
```

## Arguments

- x:

  A `VISTA` object containing DEG summaries.

- sample_comparisons:

  Optional character vector of comparison names to display.

- show_counts:

  Logical; overlay numeric counts above bars when `TRUE`.

- font_size:

  Numeric base font size for the plot.

- fill_colors:

  Named vector giving fill colors for `"Up"` and `"Down"` bars.

- facet_by:

  Either `"none"`, `"regulation"`, or `"sample_comparisons"` (alias
  `"sample_comparison"`) describing the facet variable. Use `"none"` for
  a single panel.

- facet_scale:

  Scale option passed to `facet_wrap()` when faceting.

## Value

An object returned by this function.

## Examples

``` r
NULL
#> NULL
```
