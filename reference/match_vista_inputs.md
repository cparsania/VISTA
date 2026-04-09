# Match count and metadata inputs for VISTA

`match_vista_inputs()` aligns standardized counts and sample metadata so
they can be passed directly to [`create_vista()`](create_vista.md). It
accepts the raw output from
[`read_vista_counts()`](read_vista_counts.md) or a count data
frame/matrix plus sample metadata.

## Usage

``` r
match_vista_inputs(
  counts,
  sample_info,
  column_geneid = NULL,
  sample_column = NULL,
  reorder = TRUE,
  drop_unmatched = FALSE,
  verbose = TRUE
)
```

## Arguments

- counts:

  Standardized counts from
  [`read_vista_counts()`](read_vista_counts.md), or a compatible count
  matrix/data frame.

- sample_info:

  Sample metadata from [`read_vista_metadata()`](read_vista_metadata.md)
  or a data frame coercible to that format.

- column_geneid:

  Optional gene identifier column for raw tabular counts. Ignored when
  `counts` is the list output of
  [`read_vista_counts()`](read_vista_counts.md).

- sample_column:

  Optional sample identifier column in `sample_info`.

- reorder:

  Logical; if `TRUE` (default), reorder `sample_info` to match count
  columns.

- drop_unmatched:

  Logical; if `TRUE`, keep only the intersection of count samples and
  metadata samples. Otherwise mismatches raise an error.

- verbose:

  Logical; print an informational alignment summary.

## Value

A list with standardized `counts`, aligned `sample_info`,
`column_geneid`, `sample_names`, `sample_name_map`, `row_data`, and a
small `report`.

## Examples

``` r
data("count_data", package = "VISTA")
data("sample_metadata", package = "VISTA")

cnt <- read_vista_counts(
  count_data[seq_len(25), ],
  format = "matrix",
  gene_id_column = "gene_id",
  verbose = FALSE
)
si <- read_vista_metadata(
  sample_metadata[sample_metadata$sample_names %in% cnt$sample_names, ],
  verbose = FALSE
)
matched <- match_vista_inputs(cnt, si, verbose = FALSE)

matched$column_geneid
#> [1] "gene_id"
identical(matched$sample_info$sample_names, colnames(matched$counts)[-1])
#> [1] TRUE
```
