# Read and standardize sample metadata for VISTA

`read_vista_metadata()` standardizes a sample sheet for use as
`sample_info` in
[`create_vista()`](https://cparsania.github.io/VISTA/reference/create_vista.md).
It infers or creates the required `sample_names` column using the same
conventions VISTA already accepts in the constructor.

## Usage

``` r
read_vista_metadata(
  x,
  sample_column = NULL,
  required_columns = NULL,
  drop_empty = TRUE,
  standardize_names = TRUE,
  verbose = TRUE
)
```

## Arguments

- x:

  Sample metadata as a data frame or file path.

- sample_column:

  Optional column to use as `sample_names`. If omitted, VISTA uses an
  existing `sample_names` column, non-default rownames, or common
  aliases such as `sample`, `sample_id`, or `Run`.

- required_columns:

  Optional character vector of columns that must be present after
  import.

- drop_empty:

  Logical; if `TRUE`, remove columns that are entirely `NA` or empty
  strings.

- standardize_names:

  Logical; if `TRUE`, coerce the final `sample_names` column to
  character and set rownames to match it.

- verbose:

  Logical; print an informational import summary.

## Value

A data frame suitable for use as `sample_info` in
[`create_vista()`](https://cparsania.github.io/VISTA/reference/create_vista.md).

## Examples

``` r
data("sample_metadata", package = "VISTA")

si <- read_vista_metadata(sample_metadata[seq_len(6), ])
#> Imported 6 samples with 13 metadata columns.
head(si$sample_names)
#> [1] "SRR1039508" "SRR1039509" "SRR1039512" "SRR1039513" "SRR1039516"
#> [6] "SRR1039517"
```
