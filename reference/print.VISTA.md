# Print a VISTA object like a SummarizedExperiment

Forwards to SummarizedExperiment's `show()` so the output is identical
to a plain SE. Invisibly returns `x`.

## Usage

``` r
# S3 method for class 'VISTA'
print(x, ...)
```

## Arguments

- x:

  A VISTA object.

- ...:

  Ignored.

## Value

The input object `x`, returned invisibly.

## Examples

``` r
v <- example_vista()
print(v)
#> class: SummarizedExperiment 
#> dim: 123 6 
#> metadata(12): de_results de_summary ... design comparison
#> assays(1): norm_counts
#> rownames(123): ENSG00000000003 ENSG00000000419 ... ENSG00000006607
#>   ENSG00000006625
#> rowData names(1): baseMean
#> colnames(6): SRR1039508 SRR1039512 ... SRR1039513 SRR1039517
#> colData names(14): SampleName cell ... sizeFactor sample_names
```
