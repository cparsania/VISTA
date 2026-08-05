# Barplot of DEG counts (Up/Down) across comparisons

Barplot of DEG counts (Up/Down) across comparisons

## Usage

``` r
get_deg_count_barplot(
  x,
  sample_comparisons = NULL,
  label = TRUE,
  base_size = 12,
  colors = c(Up = "red4", Down = "blue4"),
  facet_by = c("none", "regulation", "comparison"),
  facet_scales = "fixed"
)
```

## Arguments

- x:

  A `VISTA` object containing DEG summaries.

- sample_comparisons:

  Optional character vector of comparison names to display.

- label:

  Logical; overlay numeric counts above bars when `TRUE`.

- base_size:

  Numeric base font size for the plot.

- colors:

  Named vector giving fill colors for `"Up"` and `"Down"` bars.

- facet_by:

  Either `"none"`, `"regulation"`, or `"comparison"` describing the
  facet variable. Use `"none"` for a single panel.

- facet_scales:

  Scale option passed to `facet_wrap()` when faceting.

## Value

A `ggplot2` object with one bar per regulation class.

## Examples

``` r
v <- example_vista()
p <- get_deg_count_barplot(v)
print(p)
```
