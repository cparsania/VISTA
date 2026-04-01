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
v <- example_vista()
save_vista_data(v, what = "comparison", file = tempfile(fileext = ".csv"), format = "csv")
```
