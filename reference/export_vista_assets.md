# Export a complete VISTA asset bundle

Generates a standardized folder with selected VISTA plots, tabular
outputs, and a manifest describing all saved files.

## Usage

``` r
export_vista_assets(
  x,
  out_dir = "vista_assets",
  sample_comparison = NULL,
  display_id = NULL,
  include_plots = c("pca", "mds", "corr_heatmap", "deg_bar", "volcano", "ma",
    "expression_heatmap"),
  include_data = c("comparison", "norm_counts", "sample_info", "row_data", "deg_summary",
    "cutoffs"),
  plot_format = "png",
  width = 8,
  height = 6,
  heatmap_height = 10,
  units = "in",
  dpi = 300,
  top_n_labels = 50,
  heatmap_n_genes = 60,
  write_excel = FALSE,
  overwrite = TRUE
)
```

## Arguments

- x:

  A `VISTA` object.

- out_dir:

  Output directory for exported assets.

- sample_comparison:

  Optional comparison to use for comparison-specific outputs. Defaults
  to the first available comparison.

- display_id:

  Optional gene identifier column used in labeling for
  volcano/MA/heatmap plots.

- include_plots:

  Character vector of plot keys to export. Supported: `"pca"`, `"mds"`,
  `"corr_heatmap"`, `"deg_bar"`, `"deg_pie"`, `"deg_donut"`,
  `"volcano"`, `"ma"`, `"expression_heatmap"`.

- include_data:

  Character vector of data keys passed to
  [`save_vista_data()`](save_vista_data.md).

- plot_format:

  Plot format (e.g. `"png"` or `"pdf"`).

- width:

  Base plot width.

- height:

  Base plot height.

- heatmap_height:

  Height used specifically for expression heatmap export.

- units:

  Plot dimension units.

- dpi:

  Raster resolution for plots.

- top_n_labels:

  Number of top genes to annotate in MA plots.

- heatmap_n_genes:

  Number of top genes used in exported expression heatmaps.

- write_excel:

  Logical; if `TRUE`, also writes a combined XLSX workbook for all
  requested `include_data` tables (requires writexl).

- overwrite:

  Logical; if `FALSE`, aborts when `out_dir` already contains files.

## Value

Invisibly, a list with `out_dir`, `sample_comparison`, `manifest`,
`plot_files`, and `data_files`.

## Examples

``` r
v <- example_vista()
out_dir <- file.path(tempdir(), "vista_assets_example")
res <- export_vista_assets(
  v,
  out_dir = out_dir,
  include_plots = "pca",
  include_data = "comparison"
)
names(res)
#> [1] "out_dir"           "sample_comparison" "manifest"         
#> [4] "plot_files"        "data_files"       
```
