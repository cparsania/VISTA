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

  A plot object. Typically `ggplot`, `patchwork`, `Heatmap`, or
  `HeatmapList`.

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
p <- get_pca_plot(v)
out_file <- tempfile(fileext = ".pdf")
save_vista_plot(p, file = out_file, width = 7, height = 5, units = "in")
# }
```
