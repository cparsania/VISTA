# Deep validation of VISTA differential-expression fidelity

Combines structural validation from
[`validate_vista`](validate_vista.md) with backend-to-backend
equivalence checks against standalone DESeq2, edgeR, and limma runs.

## Usage

``` r
validate_vista_deep(
  counts,
  sample_info,
  column_geneid,
  group_column,
  group_numerator,
  group_denominator,
  methods = c("deseq2", "edger", "limma"),
  min_counts = 10,
  min_replicates = 1,
  log2fc_cutoff = 1,
  pval_cutoff = 0.05,
  p_value_type = "padj",
  covariates = NULL,
  design_formula = NULL,
  tolerance = 1e-08,
  return_plots = FALSE,
  error = TRUE
)
```

## Arguments

- counts:

  Raw counts with a gene-id column and sample columns.

- sample_info:

  Sample metadata.

- column_geneid:

  Column name in `counts` containing gene identifiers.

- group_column:

  Grouping column in `sample_info`.

- group_numerator:

  Numerator group(s) for pairwise comparisons.

- group_denominator:

  Denominator group(s) for pairwise comparisons.

- methods:

  Subset of `c("deseq2", "edger", "limma")` to benchmark.

- min_counts:

  Minimum total counts per gene to retain.

- min_replicates:

  Minimum samples per group meeting filtering criteria.

- log2fc_cutoff:

  Absolute log2 fold-change threshold for DEG calling.

- pval_cutoff:

  P-value or adjusted p-value threshold.

- p_value_type:

  One of `"padj"` or `"pvalue"`.

- covariates:

  Optional covariates included in the design.

- design_formula:

  Optional design formula overriding automatic design construction.

- tolerance:

  Floating-point tolerance for numerical comparisons.

- return_plots:

  When `TRUE`, include paired VISTA/reference plots for MA, volcano, DEG
  count, and PCA checks.

- error:

  When `TRUE`, abort if any discrepancy is detected.

## Value

Invisibly returns the full benchmark report from
[`benchmark_vista_equivalence`](benchmark_vista_equivalence.md).

## Examples
