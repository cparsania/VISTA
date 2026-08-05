# Pie chart of DEG counts (Up/Down) across comparisons

Pie chart of DEG counts (Up/Down) across comparisons

## Usage

``` r
get_deg_count_pieplot(
  x,
  sample_comparisons = NULL,
  label_type = c("both", "count", "percent", "none"),
  label_digits = 1,
  base_size = 12,
  colors = c(Up = "red4", Down = "blue4"),
  show_other = FALSE,
  other_color = "grey70",
  text_color = "black",
  facet_by = c("comparison", "none"),
  ncol = NULL,
  label = NULL
)
```

## Arguments

- x:

  A `VISTA` object containing DEG summaries.

- sample_comparisons:

  Optional character vector of comparison names to display.

- label_type:

  Label content: `"both"`, `"count"`, `"percent"` or `"none"`.

- label_digits:

  Integer number of decimals used for percentage labels.

- base_size:

  Numeric base font size for the plot.

- colors:

  Named vector giving fill colors for `"Up"` and `"Down"` slices. When
  `show_other = TRUE`, an `"Other"` entry may also be supplied.

- show_other:

  Logical; when `TRUE`, include non-DE genes as an `"Other"` slice using
  `other_color` unless overridden in `colors`.

- other_color:

  Fill colour used for the `"Other"` slice when `show_other = TRUE` and
  `colors` does not include `"Other"`.

- text_color:

  Colour used for pie label text.

- facet_by:

  Either `"comparison"` (default) to draw one pie per comparison, or
  `"none"` for a single pie.

- ncol:

  Optional number of columns when faceting.

- label:

  Deprecated; use `label_type`. `TRUE`/`FALSE` map to `"both"`/`"none"`.

## Value

A `ggplot2` object drawing DEG counts as a pie chart.

## Examples

``` r
v <- example_vista()
p <- get_deg_count_pieplot(v)
print(p)
```
