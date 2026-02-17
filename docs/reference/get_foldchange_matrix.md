# Extract a log2 fold-change matrix

Returns a gene-by-comparison matrix of log2 fold changes stored in a
VISTA object.

## Usage

``` r
get_foldchange_matrix(x, sample_comparisons = NULL, genes = NULL)
```

## Arguments

- x:

  A VISTA object containing differential expression results.

- sample_comparisons:

  Optional character vector of comparison names. Defaults to all
  available comparisons.

- genes:

  Optional character vector of gene identifiers. When omitted, all genes
  present in `row_data(x)` are returned.

## Value

A numeric matrix with genes in rows and comparisons in columns.
