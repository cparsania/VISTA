# Deep validation of VISTA differential-expression fidelity

This validator combines
[`validate_vista()`](https://cparsania.github.io/VISTA/reference/validate_vista.md)
with backend-to-backend numerical equivalence checks against standalone
DESeq2, edgeR, and limma runs.

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

  Raw counts (matrix/data.frame) with a gene-id column and sample
  columns.

- sample_info:

  Data frame with sample metadata.

- column_geneid:

  Column name in `counts` that contains gene identifiers.

- group_column:

  Column in `sample_info` used to group samples.

- group_numerator:

  Character vector of numerator groups for pairwise comparisons.

- group_denominator:

  Character vector of denominator groups.

- methods:

  Character vector of backends to benchmark. Any subset of
  `c("deseq2", "edger", "limma")`.

- min_counts:

  Minimum total counts per gene to retain.

- min_replicates:

  Minimum samples per group meeting filtering criteria.

- log2fc_cutoff:

  Absolute log2 fold-change threshold for DEG calling.

- pval_cutoff:

  P-value (or adjusted p-value) threshold.

- p_value_type:

  Either `"padj"` or `"pvalue"`.

- covariates:

  Optional character vector of additional sample_info columns.

- design_formula:

  Optional model formula (or formula string) including `group_column`.

- tolerance:

  Numeric tolerance used for floating-point comparisons.

- return_plots:

  Logical; if `TRUE`, return paired VISTA/reference plots for MA,
  volcano, DEG count, and PCA views.

- error:

  Logical; if `TRUE`, abort when any discrepancy is detected.

## Value

Invisibly returns the full benchmark report.

## Examples

``` r
v <- example_vista()
si <- as.data.frame(sample_info(v))
data("count_data", package = "VISTA")
count_subset <- count_data[seq_len(500), c("gene_id", si$sample_names), drop = FALSE]

report <- validate_vista_deep(
  counts = count_subset,
  sample_info = si,
  column_geneid = "gene_id",
  group_column = "cond_long",
  group_numerator = "treatment1",
  group_denominator = "control",
  methods = "limma",
  min_counts = 5,
  min_replicates = 1,
  error = FALSE
)
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes

report$valid
#> [1] TRUE

# \donttest{
data("count_data", package = "VISTA")
data("sample_metadata", package = "VISTA")

target_groups <- c("control", "treatment1")
sample_subset <- sample_metadata[sample_metadata$cond_long %in% target_groups, ]
count_subset <- count_data[seq_len(150), c("gene_id", sample_subset$sample_names)]

validate_vista_deep(
  counts = count_subset,
  sample_info = sample_subset,
  column_geneid = "gene_id",
  group_column = "cond_long",
  group_numerator = "treatment1",
  group_denominator = "control",
  methods = c("deseq2", "edger"),
  min_counts = 5,
  min_replicates = 1
)
#> estimating size factors
#> estimating dispersions
#> gene-wise dispersion estimates
#> mean-dispersion relationship
#> final dispersion estimates
#> fitting model and testing
#> estimating size factors
#> estimating dispersions
#> gene-wise dispersion estimates
#> mean-dispersion relationship
#> final dispersion estimates
#> fitting model and testing
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
# }
```
