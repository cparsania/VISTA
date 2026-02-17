# Categorize Differential Expression Results (DESeq2 or edgeR)

This function classifies genes as "Up", "Down", or "Other" based on
fold-change and p-value thresholds. It supports input from DESeq2
(`DESeqResults`), edgeR (e.g., `topTags()`), or a general tibble/data
frame.

## Usage

``` r
.categorize_deg_results(
  de_results,
  log2fc_cutoff = 1,
  pval_cutoff = 0.05,
  p_value_type = c("padj", "pvalue", "FDR", "PValue")
)
```

## Arguments

- de_results:

  A differential expression results object. Can be a DESeqResults
  object, an edgeR `topTags()` table, or a data frame/tibble containing
  fold changes and p-values.

- log2fc_cutoff:

  Numeric; minimum absolute log2 fold change to classify a DEG. Default
  is 1.

- pval_cutoff:

  Numeric; p-value or FDR threshold to define significance. Default is
  0.05.

- p_value_type:

  String; which column to use for filtering significance. One of
  `"padj"`, `"pvalue"`, `"FDR"`, `"PValue"`. Default is `"padj"`.

## Value

A tibble with gene IDs, fold changes, p-values, and a new `regulation`
column, indicating the category of each gene:

- `"Up"`: Genes with log2 fold change \>= `log2fc_cutoff` and p-value
  \<= `pval_cutoff`.

- `"Down"`: Genes with log2 fold change \<= -`log2fc_cutoff` and p-value
  \<= `pval_cutoff`.

- `"Other"`: Genes not meeting the above criteria.
