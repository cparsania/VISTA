# Coerce SummarizedExperiment to VISTA

Convert a pre-existing `SummarizedExperiment` into a valid `VISTA`
object.

## Usage

``` r
as_vista(
  se,
  assay_name = "norm_counts",
  group_column,
  comparisons = list(),
  deg_summary = list(),
  cutoffs = list(),
  group_palette = "Dark 3",
  validate = TRUE
)
```

## Arguments

- se:

  A `SummarizedExperiment` object.

- assay_name:

  Assay in `se` to use as `norm_counts` (default: `"norm_counts"`).

- group_column:

  Grouping column in `colData(se)`.

- comparisons:

  Optional named list of differential expression tables.

- deg_summary:

  Optional named list of DEG summary tables.

- cutoffs:

  Optional named list of threshold metadata.

- group_palette:

  Palette name for group colors.

- validate:

  Logical; run object validation before returning.

## Value

A `VISTA` object.

## Examples

``` r
mat <- matrix(rnorm(60), nrow = 10)
rownames(mat) <- paste0("gene", seq_len(nrow(mat)))
colnames(mat) <- paste0("sample", seq_len(ncol(mat)))
se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(norm_counts = mat),
  colData = S4Vectors::DataFrame(
    cond = rep(c("A", "B"), each = 3),
    row.names = colnames(mat)
  ),
  rowData = S4Vectors::DataFrame(
    gene_id = rownames(mat),
    row.names = rownames(mat)
  )
)
v <- as_vista(se, group_column = "cond")
v
#> class: SummarizedExperiment 
#> dim: 10 6 
#> metadata(6): de_results de_summary ... provenance vista_schema_version
#> assays(1): norm_counts
#> rownames(10): gene1 gene2 ... gene9 gene10
#> rowData names(1): gene_id
#> colnames(6): sample1 sample2 ... sample5 sample6
#> colData names(2): cond sample_names
```
