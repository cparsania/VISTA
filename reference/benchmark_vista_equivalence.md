# Benchmark VISTA against standalone differential-expression backends

Run VISTA and direct DESeq2, edgeR, and limma pipelines with matched
preprocessing, contrast definitions, and DEG thresholds, then compare
the resulting tables, DEG calls, normalized matrices, and critical plot
inputs.

## Usage

``` r
benchmark_vista_equivalence(
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
  return_plots = FALSE
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

## Value

A list with fields `valid`, `comparison_summary`, `visual_summary`, and
`methods`. Each element of `methods` contains the VISTA object, direct
backend results, structural validation output, self-consistency checks,
and optional plot objects.

## Examples

``` r
v <- example_vista()
si <- as.data.frame(sample_info(v))
data("count_data", package = "VISTA")
count_subset <- count_data[seq_len(500), c("gene_id", si$sample_names), drop = FALSE]

bm <- benchmark_vista_equivalence(
  counts = count_subset,
  sample_info = si,
  column_geneid = "gene_id",
  group_column = "cond_long",
  group_numerator = "treatment1",
  group_denominator = "control",
  methods = "limma",
  min_counts = 5,
  min_replicates = 1
)

bm$comparison_summary
#> # A tibble: 1 × 20
#>   method comparison            n_genes up_genes_identical down_genes_identical
#>   <chr>  <chr>                   <int> <lgl>              <lgl>               
#> 1 limma  treatment1_VS_control     442 TRUE               TRUE                
#> # ℹ 15 more variables: deg_sets_identical <lgl>, regulation_identical <lgl>,
#> #   norm_counts_identical <lgl>, baseMean_within_tolerance <lgl>,
#> #   log2fc_within_tolerance <lgl>, pvalue_within_tolerance <lgl>,
#> #   padj_within_tolerance <lgl>, max_abs_norm_counts_diff <dbl>,
#> #   max_abs_baseMean_diff <dbl>, max_abs_log2fc_diff <dbl>,
#> #   max_abs_pvalue_diff <dbl>, max_abs_padj_diff <dbl>, all_checks_pass <lgl>,
#> #   structural_valid <lgl>, self_consistency_valid <lgl>

# \donttest{
data("count_data", package = "VISTA")
data("sample_metadata", package = "VISTA")

target_groups <- c("control", "treatment1")
sample_subset <- sample_metadata[sample_metadata$cond_long %in% target_groups, ]
count_subset <- count_data[seq_len(150), c("gene_id", sample_subset$sample_names)]

bm <- benchmark_vista_equivalence(
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

bm$comparison_summary
#> # A tibble: 2 × 20
#>   method comparison            n_genes up_genes_identical down_genes_identical
#>   <chr>  <chr>                   <int> <lgl>              <lgl>               
#> 1 deseq2 treatment1_VS_control     123 TRUE               TRUE                
#> 2 edger  treatment1_VS_control     130 TRUE               TRUE                
#> # ℹ 15 more variables: deg_sets_identical <lgl>, regulation_identical <lgl>,
#> #   norm_counts_identical <lgl>, baseMean_within_tolerance <lgl>,
#> #   log2fc_within_tolerance <lgl>, pvalue_within_tolerance <lgl>,
#> #   padj_within_tolerance <lgl>, max_abs_norm_counts_diff <dbl>,
#> #   max_abs_baseMean_diff <dbl>, max_abs_log2fc_diff <dbl>,
#> #   max_abs_pvalue_diff <dbl>, max_abs_padj_diff <dbl>, all_checks_pass <lgl>,
#> #   structural_valid <lgl>, self_consistency_valid <lgl>
bm$visual_summary
#> # A tibble: 8 × 7
#>   method comparison  visual pass  detail structural_valid self_consistency_valid
#>   <chr>  <chr>       <chr>  <lgl> <chr>  <lgl>            <lgl>                 
#> 1 deseq2 treatment1… ma     TRUE  ident… TRUE             TRUE                  
#> 2 deseq2 treatment1… volca… TRUE  ident… TRUE             TRUE                  
#> 3 deseq2 NA          deg_c… TRUE  ident… TRUE             TRUE                  
#> 4 deseq2 NA          pca    TRUE  ident… TRUE             TRUE                  
#> 5 edger  treatment1… ma     TRUE  ident… TRUE             TRUE                  
#> 6 edger  treatment1… volca… TRUE  ident… TRUE             TRUE                  
#> 7 edger  NA          deg_c… TRUE  ident… TRUE             TRUE                  
#> 8 edger  NA          pca    TRUE  ident… TRUE             TRUE                  
# }
```
