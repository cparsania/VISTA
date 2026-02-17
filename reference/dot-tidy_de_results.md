# Standardize Differential Expression Results

Harmonize DE results (DESeq2/edgeR) to a common schema and ensure gene
IDs are stored as rownames. Returns S4Vectors::DataFrame.

## Usage

``` r
.tidy_de_results(tbl, rowname_col = NULL)
```

## Arguments

- tbl:

  data.frame/tibble or S4Vectors::DataFrame with DE results.

- rowname_col:

  optional column to promote to rownames if missing.

## Value

S4Vectors::DataFrame with rownames = gene IDs and columns like log2FC,
pvalue, padj (if present), baseMean/AveExpr mapped to baseMean, etc.
