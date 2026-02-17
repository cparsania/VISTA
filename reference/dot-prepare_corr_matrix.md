# Prepare correlation matrix from normalized expression

Prepare correlation matrix from normalized expression

## Usage

``` r
.prepare_corr_matrix(mat, meta, corr_method = "pearson")
```

## Arguments

- mat:

  A numeric matrix of log-normalized gene expression (genes × samples).

- corr_method:

  Correlation method; one of "pearson", "kendall", or "spearman".

## Value

A correlation matrix.
