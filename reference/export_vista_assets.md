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
# \donttest{
set.seed(1)
mat <- matrix(rpois(60, lambda = 20), nrow = 10)
rownames(mat) <- paste0("gene", seq_len(nrow(mat)))
colnames(mat) <- paste0("sample", seq_len(ncol(mat)))
se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(norm_counts = mat),
  colData = S4Vectors::DataFrame(
    group = rep(c("ctrl", "trt"), each = 3),
    row.names = colnames(mat)
  ),
  rowData = S4Vectors::DataFrame(
    gene_id = rownames(mat),
    row.names = rownames(mat)
  )
)
de <- data.frame(
  gene_id = rownames(mat),
  log2fc = rnorm(nrow(mat)),
  pvalue = runif(nrow(mat)),
  padj = runif(nrow(mat)),
  regulation = "Other",
  row.names = rownames(mat)
)
v <- as_vista(se, group_column = "group")
md <- S4Vectors::metadata(v)
md$de_results <- S4Vectors::SimpleList(trt_vs_ctrl = de)
md$de_summary <- S4Vectors::SimpleList(trt_vs_ctrl = as.data.frame(table(de$regulation)))
S4Vectors::metadata(v) <- md
out_dir <- tempfile("vista_assets_")
export_vista_assets(
  v,
  out_dir = out_dir,
  include_plots = "pca",
  include_data = c("comparison", "norm_counts")
)
# }
```
