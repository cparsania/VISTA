# Extract a log2 fold-change matrix

Returns a gene-by-comparison matrix of log2 fold changes stored in a
VISTA object.

## Usage

``` r
get_foldchange_matrix(
  x,
  sample_comparisons = NULL,
  genes = NULL,
  display_id = NULL,
  display_from = NULL,
  display_orgdb = NULL
)
```

## Arguments

- x:

  A VISTA object containing differential expression results.

- sample_comparisons:

  Optional character vector of comparison names. Defaults to all
  available comparisons.

- genes:

  Optional character vector of gene identifiers. When omitted, all genes
  present in `row_data(x)` are returned. May be given as display labels
  when `display_id` is set.

- display_id:

  Optional `rowData()` column (or annotation key) naming the labels to
  work in, matching
  [`get_foldchange_heatmap()`](https://cparsania.github.io/VISTA/reference/get_foldchange_heatmap.md).
  When supplied, `genes` is accepted in those labels and the returned
  matrix is labelled with them. Identifiers with no label keep their
  original value.

- display_from:

  Identifier type `rownames(x)` are in, used only when `display_id` is
  not a `rowData()` column and the mapping goes through an annotation
  package.

- display_orgdb:

  An `OrgDb` used for that mapping.

## Value

A numeric matrix with genes in rows and comparisons in columns. Rows
keep the object's gene identifiers unless `display_id` is supplied.
Because display labels need not be unique, duplicates are made unique
with a warning so that every row stays addressable by name.

## Examples

``` r
v <- example_vista()
mat <- get_foldchange_matrix(v)
dim(mat)
#> [1] 123   1
```
