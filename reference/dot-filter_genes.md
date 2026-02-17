# Filter genes by user-specified IDs or variability

This internal function subsets the normalized expression matrix to
retain only selected genes (by name) or top variable genes (by
variance).

## Usage

``` r
.filter_genes(mat, genes = NULL, top_n_genes = NULL)
```

## Arguments

- mat:

  A numeric matrix with genes as rows and samples as columns.

- genes:

  Optional character vector of gene IDs to retain.

- top_n_genes:

  Optional integer; the number of top variable genes to keep.

## Value

A filtered matrix.
