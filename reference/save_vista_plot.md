# Save a VISTA plot object to disk

Saves plot objects returned by VISTA plotting functions to file.
Supports both `ggplot`-like objects (saved via
[`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html))
and `ComplexHeatmap` objects (`Heatmap` / `HeatmapList`) saved via
graphics devices.

## Usage

``` r
save_vista_plot(
  plot,
  file,
  width = 8,
  height = 6,
  units = "in",
  dpi = 300,
  device = NULL,
  ...
)
```

## Arguments

- plot:

  A plot object. Typically `ggplot`, `patchwork`, `Heatmap`,
  `HeatmapList`, or a `recordedplot` (as returned by
  `get_enrichment_chord(return_type = "plot")`).

- file:

  Output file path.

- width:

  Plot width.

- height:

  Plot height.

- units:

  Units for `width` and `height`. One of `"in"`, `"cm"`, `"mm"`, or
  `"px"`.

- dpi:

  Resolution for raster outputs.

- device:

  Optional graphics device (e.g. `"png"`, `"pdf"`). If `NULL`, inferred
  from `file` extension (defaults to `"png"` when missing).

- ...:

  Additional arguments passed to
  [`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html)
  for ggplot-like objects.

## Value

Invisibly, the normalized output file path.

## Examples

``` r
v <- example_vista()
p <- get_pca_plot(v)
out_file <- tempfile(fileext = ".pdf")
save_vista_plot(p, file = out_file, width = 7, height = 5, units = "in")
```
