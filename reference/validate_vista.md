# Validate a VISTA object

Performs structural and metadata contract checks on a `VISTA` object.
Use `level = "core"` for invariants that all VISTA objects must satisfy,
or `level = "full"` to additionally validate method-specific metadata
used by advanced workflows.

## Usage

``` r
validate_vista(x, level = c("core", "full"), error = TRUE)
```

## Arguments

- x:

  A `VISTA` object.

- level:

  Validation depth. One of `"core"` or `"full"`.

- error:

  Logical; if `TRUE`, aborts on validation failures. If `FALSE`, returns
  the validation report and emits a warning when issues are found.

## Value

Invisibly returns a list with fields `valid`, `level`, and `issues`.

## Examples

``` r
mat <- matrix(rnorm(60), nrow = 10)
rownames(mat) <- paste0("gene", seq_len(nrow(mat)))
colnames(mat) <- paste0("sample", seq_len(ncol(mat)))
se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(norm_counts = mat),
  colData = S4Vectors::DataFrame(
    cond = rep(c("A", "B"), each = 3),
    row.names = colnames(mat)
  ),
  rowData = S4Vectors::DataFrame(
    gene_id = rownames(mat),
    row.names = rownames(mat)
  )
)
v <- as_vista(se, group_column = "cond")
validate_vista(v)
```
