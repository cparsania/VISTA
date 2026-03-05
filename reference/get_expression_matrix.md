# Retrieve an expression matrix from a VISTA object

Returns the specified assay (default `"norm_counts"`) with optional
gene/sample subsetting, replicate summarization by group, and simple
transformations.

## Usage

``` r
get_expression_matrix(
  x,
  assay_name = "norm_counts",
  genes = NULL,
  sample_names = NULL,
  group_column = NULL,
  summarise = FALSE,
  transform = c("none", "log2", "zscore")
)
```

## Arguments

- x:

  A `VISTA` object.

- assay_name:

  Assay to extract (default: `"norm_counts"`).

- genes:

  Optional character vector of gene IDs to keep.

- sample_names:

  Optional character vector of sample IDs (or group labels when
  `summarise = TRUE`) to keep.

- group_column:

  Optional column in `sample_info` used when summarising replicates.
  Defaults to the stored `group_column`.

- summarise:

  Logical; if `TRUE`, averages replicates per `group_column`.

- transform:

  One of `"none"`, `"log2"`, or `"zscore"` applied after subsetting and
  optional summarisation. Z-scores are computed within the selected
  samples/columns.

## Value

A numeric matrix with genes in rows and samples (or groups) in columns.

## Examples

``` r
NULL
#> NULL
```
