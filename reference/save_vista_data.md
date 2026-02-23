# Save VISTA tabular outputs to disk

Exports selected data components from a `VISTA` object to
CSV/TSV/RDS/XLSX.

## Usage

``` r
save_vista_data(
  x,
  what = c("comparison", "comparisons", "norm_counts", "sample_info", "row_data",
    "deg_summary", "cutoffs"),
  file,
  sample_comparison = NULL,
  format = NULL,
  include_rownames = TRUE
)
```

## Arguments

- x:

  A `VISTA` object.

- what:

  Character vector specifying which object(s) to export. Supported
  values are `"comparison"`, `"comparisons"`, `"norm_counts"`,
  `"sample_info"`, `"row_data"`, `"deg_summary"`, and `"cutoffs"`.

- file:

  Output file path.

- sample_comparison:

  Optional comparison name used when `what` includes `"comparison"`.
  Defaults to the first comparison in `comparisons(x)`.

- format:

  Output format. One of `"csv"`, `"tsv"`, `"rds"`, `"xlsx"`. If `NULL`,
  inferred from `file` extension.

- include_rownames:

  Logical; include meaningful row identifiers (e.g., gene IDs or sample
  names) as explicit columns where applicable.

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
save_vista_data(v, what = "comparison", file = tempfile(fileext = ".csv"), format = "csv")
save_vista_data(
  v,
  what = c("comparison", "norm_counts"),
  file = tempfile(fileext = ".rds"),
  format = "rds"
)
# }
```
