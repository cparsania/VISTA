# Build a small example VISTA object

Creates a lightweight `VISTA` object from built-in package datasets
(`count_data` and `sample_metadata`) for use in examples and tutorials.
The default call returns a precomputed object to keep examples and
package checks fast. Non-default argument combinations fall back to
rebuilding the object from the packaged example inputs.

## Usage

``` r
example_vista(n_genes = 150, n_per_group = 3, method = "deseq2")
```

## Arguments

- n_genes:

  Number of genes to include (default `150`).

- n_per_group:

  Number of samples per group (`control` and `treatment1`) to include
  (default `3`).

- method:

  Differential expression backend passed to
  [`create_vista()`](create_vista.md).

## Value

A `VISTA` object.

## Examples

``` r
v <- example_vista()
v
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
