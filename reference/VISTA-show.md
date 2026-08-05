# Display a VISTA object

Prints the `SummarizedExperiment` header followed by a short summary of
the analysis state – grouping, comparisons, active
differential-expression source, thresholds and schema version – none of
which was visible before.

## Usage

``` r
# S4 method for class 'VISTA'
show(object)
```

## Arguments

- object:

  A `VISTA` object.

## Value

`object`, invisibly. Called for its output.

## Examples

``` r
v <- example_vista()
v
#> class: VISTA 
#> dim: 123 6 
#> metadata(12): de_results de_summary ... design comparison
#> assays(2): norm_counts counts
#> rownames(123): ENSG00000000003 ENSG00000000419 ... ENSG00000006607
#>   ENSG00000006625
#> rowData names(1): baseMean
#> colnames(6): SRR1039508 SRR1039512 ... SRR1039513 SRR1039517
#> colData names(14): SampleName cell ... sizeFactor sample_names
#> -------- VISTA --------
#> group column: cond_long (control, treatment1)
#> comparisons: treatment1_VS_control
#> DE source: deseq2
#> cutoffs: |log2FC| >= 1, padj <= 0.05
#> raw counts: available via counts()
#> schema: 1.1.0
```
