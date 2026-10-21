# Run Differential Expression Analysis with DESeq2, edgeR, or limma-voom

These functions encapsulate the standard RNA-seq analysis workflow using
DESeq2 (run_deseq_analysis), edgeR (run_edger_analysis), or limma-voom
(run_limma_analysis), including: gene filtering, design matrix setup,
normalization, model fitting, differential testing, DEG classification
(`"Up"`, `"Down"`, `"Other"`), and result formatting.

Both methods return output in a harmonized structure ready for
downstream use in
[create_vista](https://cparsania.github.io/VISTA/reference/create_vista.md)
or standalone DEG summaries.

## Usage

``` r
run_deseq_analysis(
  counts,
  sample_info,
  column_geneid,
  group_column,
  group_numerator,
  group_denominator,
  covariates = NULL,
  design_formula = NULL,
  min_counts = 10,
  min_replicates = 1,
  log2fc_cutoff = 1,
  pval_cutoff = 0.05,
  p_value_type = "padj"
)

run_edger_analysis(
  counts,
  sample_info,
  column_geneid,
  group_column,
  group_numerator,
  group_denominator,
  covariates = NULL,
  design_formula = NULL,
  min_counts = 10,
  min_replicates = 1,
  log2fc_cutoff = 1,
  pval_cutoff = 0.05,
  p_value_type = "FDR",
  edger_test = c("LRT", "QL")
)

run_limma_analysis(
  counts,
  sample_info,
  column_geneid,
  group_column,
  group_numerator,
  group_denominator,
  covariates = NULL,
  design_formula = NULL,
  min_counts = 10,
  min_replicates = 1,
  log2fc_cutoff = 1,
  pval_cutoff = 0.05,
  p_value_type = "FDR"
)
```

## Arguments

- counts:

  A data frame or matrix of raw counts with one gene per row. Must
  include a column defined by `column_geneid`, and column names must
  match entries in `sample_info$sample_names`.

- sample_info:

  A data frame with sample metadata. Must contain `sample_names` and the
  specified grouping column.

- column_geneid:

  A string identifying the column name containing gene identifiers.

- group_column:

  The name of the column in `sample_info` that defines experimental
  groups.

- group_numerator:

  A character vector of numerator group(s) for fold-change comparisons.

- group_denominator:

  A character vector of denominator group(s) for fold-change
  comparisons.

- covariates:

  Optional character vector of additional sample_info columns to adjust
  for.

- design_formula:

  Optional model formula (or formula string). When provided, it
  overrides automatic design construction from `group_column` +
  `covariates`. Must include `group_column`.

- min_counts:

  Minimum per-sample read count a gene must reach to count as detected
  in that sample. Also applied as a minimum row total across all samples
  in an initial pre-filter. Default: `10`.

- min_replicates:

  Minimum number of samples in which a gene must reach `min_counts` for
  that gene to enter the model. Counted across the whole experiment, not
  within each group. Default: `1`.

- log2fc_cutoff:

  Absolute log2 fold-change threshold to define DEGs. Default: `1`.

- pval_cutoff:

  P-value or adjusted p-value cutoff for significance. Default: `0.05`.

- p_value_type:

  For DESeq2: one of `"padj"` or `"pvalue"`. For edgeR/limma: one of
  `"FDR"` or `"PValue"`.

- edger_test:

  edgeR only: `"LRT"` (likelihood-ratio test,
  [glmFit](https://rdrr.io/pkg/edgeR/man/glmfit.html) and
  [glmLRT](https://rdrr.io/pkg/edgeR/man/glmLRT.html); the default) or
  `"QL"` (quasi-likelihood F-test,
  [glmQLFit](https://rdrr.io/pkg/edgeR/man/glmQLFit.html) and
  [glmQLFTest](https://rdrr.io/pkg/edgeR/man/glmQLFTest.html), as the
  edgeR authors recommend for bulk RNA-seq). The default changes to
  `"QL"` in VISTA 1.4.0; until then, leaving it unset (or `NULL`) prints
  a notice. Set it explicitly to keep results stable.

## Value

A named list with components:

- `norm_counts`: Matrix of normalized expression values (CPM for
  edgeR/limma, DESeq2-normalized counts).

- `sample_info`: Updated sample metadata.

- `row_data`: Gene-level metadata, including mean expression.

- `comparisons`: Named list of DEG result tibbles (one per comparison),
  each containing standardized columns: `gene_id`, `log2fc`, `pvalue`,
  `p.adj`, and `regulation`.

- `deg_summary`: List of summary tables showing DEG regulation counts.

## Details

Perform differential expression (DE) analysis across multiple group
comparisons using DESeq2, edgeR, or limma-voom. These functions process
raw count data, normalize it, execute pairwise group-level tests, and
return standardized DEG outputs compatible with `VISTA`-based
visualization and analysis.

- For DESeq2, normalization is performed via
  [DESeq](https://rdrr.io/pkg/DESeq2/man/DESeq.html), and DE testing
  uses [results](https://rdrr.io/pkg/DESeq2/man/results.html).

- For edgeR, normalization uses
  [calcNormFactors](https://rdrr.io/pkg/edgeR/man/calcNormFactors.html),
  and testing uses [glmLRT](https://rdrr.io/pkg/edgeR/man/glmLRT.html)
  or, with `edger_test = "QL"`,
  [glmQLFTest](https://rdrr.io/pkg/edgeR/man/glmQLFTest.html).

- For limma, normalization uses
  [calcNormFactors](https://rdrr.io/pkg/edgeR/man/calcNormFactors.html) +
  [voom](https://rdrr.io/pkg/limma/man/voom.html), and testing uses
  [eBayes](https://rdrr.io/pkg/limma/man/ebayes.html).

Low-abundance filtering is applied before model fitting, and for the
edgeR/limma backends before
[calcNormFactors](https://rdrr.io/pkg/edgeR/man/calcNormFactors.html) so
that normalization factors are not estimated from genes that will be
discarded. All three backends use the same filtering predicate: a gene
is retained when at least `min_replicates` samples reach `min_counts`.
Filtering is applied to the counts as given; DESeq2 then receives
non-integer (estimated) counts rounded to the nearest integer, as
`DESeqDataSetFromTximport()` does, with a message. Every backend refuses
negative values and warns about non-integer values with an unmistakable
sign of not being counts (every sample summing to exactly one million,
or small, almost all fractional values). VISTA cannot otherwise verify
that its input is counts: supply raw or estimated counts, not normalised
values. Gene regulation status is determined via
[`.categorize_deg_results()`](https://cparsania.github.io/VISTA/reference/dot-categorize_deg_results.md)
based on user thresholds.

All output comparison results are internally standardized via
[`.tidy_de_results()`](https://cparsania.github.io/VISTA/reference/dot-tidy_de_results.md)
to ensure a uniform column schema compatible with VISTA plotting tools.

## See also

[create_vista](https://cparsania.github.io/VISTA/reference/create_vista.md),
[DESeq](https://rdrr.io/pkg/DESeq2/man/DESeq.html),
[glmLRT](https://rdrr.io/pkg/edgeR/man/glmLRT.html),
[voom](https://rdrr.io/pkg/limma/man/voom.html)

## Examples

``` r
v <- example_vista()
si <- as.data.frame(sample_info(v))
data("count_data", package = "VISTA")
counts_small <- count_data[seq_len(200), c("gene_id", si$sample_names), drop = FALSE]
limma_results <- run_limma_analysis(
  counts = counts_small,
  sample_info = si,
  column_geneid = "gene_id",
  group_column = "cond_long",
  group_numerator = "treatment1",
  group_denominator = "control",
  min_counts = 5,
  min_replicates = 1
)
#> calcNormFactors has been renamed to normLibSizes
names(limma_results$comparisons)
#> [1] "treatment1_VS_control"
```
