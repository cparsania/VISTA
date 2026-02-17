# Perform pairwise DE comparisons using DESeq2 results

Perform pairwise DE comparisons using DESeq2 results

## Usage

``` r
.run_deseq_comparisons(
  dds,
  group_column,
  group_numerator,
  group_denominator,
  log2fc_cutoff,
  pval_cutoff,
  p_value_type,
  column_geneid
)
```

## Arguments

- dds:

  A DESeq2 object after running DESeq().

- group_column:

  The column used for grouping in the design formula.

- group_numerator:

  Vector of numerator group names.

- group_denominator:

  Vector of denominator group names.

- log2fc_cutoff:

  Absolute log2 fold-change cutoff for DEG classification.

- pval_cutoff:

  P-value threshold for DEG classification.

- p_value_type:

  Column to use for p-value filtering ("padj" or "pvalue").

- column_geneid:

  Name of the gene ID column to include in output.

## Value

A list with:

- `comparisons`: List of DEG result tibbles (with gene ID and fold
  changes)

- `deg_summary`: Summary table of DEG counts by regulation
