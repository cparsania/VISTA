# Cluster genes by their log2FC profiles using k-means

Cluster genes by their log2FC profiles using k-means

## Usage

``` r
.cluster_log2fc_matrix(df, gene_id_col = "display_gene", k = 3)
```

## Arguments

- df:

  A long data.frame with columns: gene_id_col, comparison, log2fc.

- gene_id_col:

  Character string of column name for gene IDs (e.g. "gene_name" or
  "display_gene").

- k:

  Number of clusters.

## Value

A data.frame with gene IDs and cluster assignments.
