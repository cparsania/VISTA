# Plot cell-type deconvolution heatmap

Plot cell-type deconvolution heatmap

## Usage

``` r
get_celltype_heatmap(
  x,
  group_column = NULL,
  samples = NULL,
  cell_types = NULL,
  top_n = 20,
  transform = c("none", "zscore", "log1p"),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_values = FALSE,
  font_size = 11,
  return_type = c("plot", "matrix", "both")
)
```

## Arguments

- x:

  A VISTA object.

- group_column:

  Optional grouping column from `sample_info(x)`. If provided and
  `cluster_columns = FALSE`, samples are ordered by this group.

- samples:

  Optional character vector of samples to include.

- cell_types:

  Optional character vector of cell types to include.

- top_n:

  Number of top cell types by mean score when `cell_types` is `NULL`.

- transform:

  One of `"none"`, `"zscore"`, or `"log1p"`.

- cluster_rows:

  Logical; hierarchical cluster cell types.

- cluster_columns:

  Logical; hierarchical cluster samples.

- show_values:

  Logical; overlay numeric values on tiles.

- font_size:

  Base font size.

- return_type:

  One of `"plot"`, `"matrix"`, or `"both"`.

## Value

A ggplot object, matrix, or list depending on `return_type`.

## Examples

``` r
mat <- matrix(rpois(20, lambda = 20), nrow = 5)
rownames(mat) <- paste0("gene", seq_len(5))
colnames(mat) <- paste0("sample", seq_len(4))
se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(norm_counts = mat),
  colData = S4Vectors::DataFrame(
    cond = c("A", "A", "B", "B"),
    row.names = colnames(mat)
  ),
  rowData = S4Vectors::DataFrame(
    gene_id = rownames(mat),
    row.names = rownames(mat)
  )
)
v <- as_vista(se, group_column = "cond")
md <- S4Vectors::metadata(v)
md$cell_fractions <- data.frame(
  fibroblast = c(0.2, 0.3, 0.4, 0.5),
  epithelial = c(0.8, 0.7, 0.6, 0.5),
  row.names = colnames(mat)
)
S4Vectors::metadata(v) <- md
get_celltype_heatmap(v, group_column = "cond")
```
