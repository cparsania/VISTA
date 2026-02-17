# Align a DE table to reference gene order; add NA rows for missing genes. Keeps a real `gene_id` column and returns a base data.frame. Canonical column order: gene_id, baseMean, log2fc, lfcSE, stat, pvalue, padj, regulation, ...

Align a DE table to reference gene order; add NA rows for missing genes.
Keeps a real `gene_id` column and returns a base data.frame. Canonical
column order: gene_id, baseMean, log2fc, lfcSE, stat, pvalue, padj,
regulation, ...

## Usage

``` r
.align_de_to_counts(
  df,
  ref_rn,
  id_col = NULL,
  strict = FALSE,
  warn_missing = TRUE
)
```
