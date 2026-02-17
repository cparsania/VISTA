# VISTA constructor

Create a VISTA object that extends `SummarizedExperiment`.

## Usage

``` r
vista(
  norm_counts,
  sample_info,
  row_data,
  comparisons = list(),
  deg_summary = list(),
  cutoffs = list(),
  group_column = character(1),
  group_palette = "Dark 3"
)
```

## Arguments

- norm_counts:

  A normalized expression matrix (or `data.frame`) with genes as rows
  and samples as columns. Must have rownames.

- sample_info:

  A `data.frame` with sample metadata. Row names must match
  `colnames(norm_counts)`.

- row_data:

  A `data.frame` with feature metadata. Row names must match
  `rownames(norm_counts)`.

- comparisons:

  A named list of differential expression result tables (e.g.
  [`S4Vectors::DataFrame`](https://rdrr.io/pkg/S4Vectors/man/DataFrame-class.html)s)
  keyed by comparison name. Default: empty list.

- deg_summary:

  A named list of per-comparison DEG summary tables. Default: empty
  list.

- cutoffs:

  A named list of thresholds used for DEG classification (e.g. `alpha`,
  `lfc`). Default: empty list.

- group_column:

  A **non-empty** character scalar naming the grouping column in
  `sample_info`. This is **required** and must exist in `sample_info`.

- group_palette:

  A colorspace qualitative palette name used to assign colors to groups
  in `group_column`. One of
  `c("Pastel 1","Dark 2","Dark 3","Set 2","Set 3","Warm","Cold","Harmonic","Dynamic")`.
  Default: `"Dark 3"`.

## Value

A `VISTA` object. The following keys are stored in `metadata(vista)`:

- `$de_results`: named `SimpleList` of DE result tables.

- `$de_summary`: named `SimpleList` of DEG summary tables.

- `$de_cutoffs`: named list of thresholds.

- `$group`: list with `column`, `palette`, `colors`.

- `$provenance`: list with constructor version, timestamp, session info.

## Examples

``` r
# Low-level constructor (advanced usage)
# Most users should use create_vista() instead

# Create sample data
norm_counts_mat <- matrix(
  rnorm(100), nrow = 10, ncol = 10,
  dimnames = list(
    paste0("gene", 1:10),
    paste0("sample", 1:10)
  )
)

sample_info_df <- data.frame(
  groups = rep(c("A", "B"), each = 5),
  row.names = colnames(norm_counts_mat)
)

row_data_df <- data.frame(
  gene_id = rownames(norm_counts_mat),
  row.names = rownames(norm_counts_mat)
)

# Create VISTA object
v <- vista(
  norm_counts = norm_counts_mat,
  sample_info = sample_info_df,
  row_data = row_data_df,
  group_column = "groups"
)

# Examine object
v
#> class: SummarizedExperiment 
#> dim: 10 10 
#> metadata(5): de_results de_summary de_cutoffs group provenance
#> assays(1): norm_counts
#> rownames(10): gene1 gene2 ... gene9 gene10
#> rowData names(1): gene_id
#> colnames(10): sample1 sample2 ... sample9 sample10
#> colData names(1): groups
norm_counts(v)[1:5, 1:5]
#>          sample1    sample2    sample3     sample4     sample5
#> gene1  1.7150650  1.7869131 -1.6866933  0.68864025 -1.12310858
#> gene2  0.4609162  0.4978505  0.8377870  0.55391765 -0.40288484
#> gene3 -1.2650612 -1.9666172  0.1533731 -0.06191171 -0.46665535
#> gene4 -0.6868529  0.7013559 -1.1381369 -0.30596266  0.77996512
#> gene5 -0.4456620 -0.4727914  1.2538149 -0.38047100 -0.08336907
```
