# Subset a VISTA object

Subsets the underlying `SummarizedExperiment` and keeps the analysis
metadata consistent with the result.

## Usage

``` r
# S4 method for class 'VISTA,ANY,ANY,ANY'
x[i, j, ..., drop = FALSE]
```

## Arguments

- x:

  A `VISTA` object.

- i, j:

  Row (gene) and column (sample) subscripts.

- ...:

  Passed to the `SummarizedExperiment` method.

- drop:

  Ignored, as for `SummarizedExperiment`.

## Value

A `VISTA` object.

## Details

Row (gene) subsetting reindexes every stored differential-expression
table – the active `de_results`/`de_summary` and every entry of
`de_results_by_method`/`de_summary_by_method` – to the retained genes,
in the new row order. DEG summaries are recounted from the retained
rows.

Column (sample) subsetting leaves the DE tables untouched, because they
are analysis *results* rather than per-sample data; recomputing them
would require re-running the model. Group and comparison colour maps are
pruned to the groups that still have samples, and a warning names any
comparison whose groups are no longer represented, since plots of that
comparison no longer correspond to the samples in the object.

Each subset appends an entry to `metadata(x)$provenance$subset_history`.

## Examples

``` r
v <- example_vista()

# DE tables follow the genes.
v10 <- v[seq_len(10), ]
nrow(comparisons(v10)[[1]])
#> [1] 10

# Sample subsetting prunes the group colour map.
v_ctrl <- v[, sample_info(v)$cond_long == "control"]
#> Warning: Comparison "treatment1_VS_control" is no longer represented by the retained
#> samples.
#> ℹ The stored results are kept unchanged; re-run `create_vista()` to analyse
#>   this subset.
group_colors(v_ctrl)
#>   control 
#> "#C87A8A" 
```
