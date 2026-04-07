# Create a VISTA Object with Internal DE Analysis

This wrapper performs differential expression (DE) analysis (DESeq2,
edgeR, limma, or both) and returns a fully initialized `VISTA` object.
The object stores expression matrices and annotations in the
SummarizedExperiment core, while all DE outputs and configuration live
in `metadata(vista)`:

- `$de_results`: named `SimpleList` of per-contrast DE tables

- `$de_summary`: named `SimpleList` of summary tables

- `$de_cutoffs`: list of thresholds/method options

- `$group`: list with `column`, `palette`, `colors`

## Usage

``` r
create_vista(
  counts,
  sample_info,
  column_geneid,
  group_column,
  group_numerator,
  group_denominator,
  method = c("deseq2", "edger", "limma", "both"),
  min_counts = 10,
  min_replicates = 1,
  log2fc_cutoff = 1,
  pval_cutoff = 0.05,
  p_value_type = "padj",
  covariates = NULL,
  design_formula = NULL,
  consensus_mode = c("intersection", "union"),
  consensus_log2fc = c("mean", "deseq2", "edger"),
  result_source = NULL,
  group_palette = "Dark 2",
  comparison_palette = "Dark 3",
  validate = TRUE
)
```

## Arguments

- counts:

  Raw counts (matrix/data.frame) with a gene-id column and sample
  columns.

- sample_info:

  Data frame with sample metadata. Must contain `sample_names` (or have
  rownames equal to sample columns in `counts`) and the `group_column`.

- column_geneid:

  Column name in `counts` that contains gene identifiers.

- group_column:

  Column in `sample_info` used to group samples.

- group_numerator:

  Character vector of numerator groups for pairwise comparisons.

- group_denominator:

  Character vector of denominator groups (same length/order as
  numerator).

- method:

  `"deseq2"`, `"edger"`, `"limma"`, or `"both"`.

- min_counts:

  Minimum total counts per gene to retain (default: 10).

- min_replicates:

  Minimum samples per group meeting `min_counts` (default: 1).

- log2fc_cutoff:

  Absolute LFC threshold for DEG calling (default: 1).

- pval_cutoff:

  p-value (raw or adjusted) threshold (default: 0.05).

- p_value_type:

  Which p-value column to use (`"padj"` or `"pvalue"`). Default:
  `"padj"`.

- covariates:

  Optional character vector of additional sample_info columns to adjust
  for. These are included as additive terms in the DE design.

- design_formula:

  Optional model formula (or formula string) overriding automatic
  construction from `group_column` + `covariates`. Must include
  `group_column`.

- consensus_mode:

  When `method = "both"`, how to define consensus calls:
  `"intersection"` (both methods significant in same direction) or
  `"union"` (either method significant; discordant directions excluded).

- consensus_log2fc:

  When `method = "both"`, how to populate consensus `log2fc`: `"mean"`,
  `"deseq2"`, or `"edger"`.

- result_source:

  Active DE source used in `metadata(v)$de_results`. For
  `method = "both"`, one of `"consensus"`, `"deseq2"`, `"edger"`. For
  single-method runs, this must match `method`.

- group_palette:

  Qualitative palette name for
  [`colorspace::qualitative_hcl()`](https://colorspace.R-Forge.R-project.org/reference/hcl_palettes.html).
  One of
  `c("Pastel 1","Dark 2","Dark 3","Set 2","Set 3","Warm","Cold","Harmonic","Dynamic")`.
  Default: `"Dark 2"`.

- comparison_palette:

  Qualitative palette name used to assign colors per comparison (stored
  in `metadata(v)$comparison$colors`). Defaults to `"Dark 3"`.

- validate:

  Logical; if `TRUE` (default), run full
  [`validate_vista()`](validate_vista.md) checks before returning the
  object.

## Value

A `VISTA` object:

- **assays(v)**: `norm_counts` (matrix)

- **colData(v)**: `sample_info` (DataFrame)

- **rowData(v)**: `row_data` (DataFrame)

- **metadata(v)**: `de_results`, `de_summary`, `de_cutoffs`, `group`,
  `comparison`, `provenance`

## Details

Contrast names follow `"numerator_VS_denominator"`. Each DE table must
have rownames identical to the final `norm_counts` rownames. When
`method = "both"`, method-specific and consensus DE tables are stored in
`metadata(v)$de_results_by_method` and
`metadata(v)$de_summary_by_method`, and the active source is tracked in
`metadata(v)$de_active_source`.

## See also

[as_vista](as_vista.md), [VISTA-class](VISTA-class.md),
[qualitative_hcl](https://colorspace.R-Forge.R-project.org/reference/hcl_palettes.html)

## Examples

``` r
# Load example data
data("count_data", package = "VISTA")
data("sample_metadata", package = "VISTA")

# Create VISTA object with DESeq2 (default method)
vista <- create_vista(
  counts = count_data[seq_len(100), ],
  sample_info = sample_metadata[seq_len(6), ],
  column_geneid = "gene_id",
  group_column = "cond_long",
  group_numerator = "treatment1",
  group_denominator = "control",
  log2fc_cutoff = 0.6,
  pval_cutoff = 0.05
)
#> estimating size factors
#> estimating dispersions
#> gene-wise dispersion estimates
#> mean-dispersion relationship
#> final dispersion estimates
#> fitting model and testing

# Examine the VISTA object
vista
#> class: SummarizedExperiment 
#> dim: 85 6 
#> metadata(12): de_results de_summary ... design comparison
#> assays(1): norm_counts
#> rownames(85): ENSG00000000003 ENSG00000000419 ... ENSG00000005469
#>   ENSG00000005471
#> rowData names(1): baseMean
#> colnames(6): SRR1039508 SRR1039509 ... SRR1039516 SRR1039517
#> colData names(14): SampleName cell ... sizeFactor sample_names

# Access comparisons
names(comparisons(vista))
#> [1] "treatment1_VS_control"

# View DEG summary
deg_summary(vista)
#> $treatment1_VS_control
#>   regulation  n
#> 1       Down  1
#> 2      Other 82
#> 3         Up  2
#> 

# View cutoffs used
cutoffs(vista)
#> $log2fc
#> [1] 0.6
#> 
#> $pval
#> [1] 0.05
#> 
#> $p_value_type
#> [1] "padj"
#> 
#> $method
#> [1] "deseq2"
#> 
#> $min_counts
#> [1] 10
#> 
#> $min_replicates
#> [1] 1
#> 
#> $covariates
#> character(0)
#> 
#> $design_formula
#> NULL
#> 
#> $consensus_mode
#> NULL
#> 
#> $consensus_log2fc
#> NULL
#> 
#> $active_source
#> [1] "deseq2"
#> 

# Multiple comparisons example
# \donttest{
vista_multi <- create_vista(
  counts = count_data,
  sample_info = sample_metadata,
  column_geneid = "gene_id",
  group_column = "cell",
  group_numerator = c("N052611", "N080611"),
  group_denominator = c("N61311", "N61311"),
  method = "edger",
  log2fc_cutoff = 1.0,
  pval_cutoff = 0.01
)
# }
```
