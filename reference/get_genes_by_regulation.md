# Get Genes by Regulation

Extract gene IDs by regulation class from selected comparisons in a
VISTA object.

## Usage

``` r
get_genes_by_regulation(x, sample_comparisons, regulation = "Both")
```

## Arguments

- x:

  A VISTA object.

- sample_comparisons:

  Character vector of comparison names to include.

- regulation:

  One of "Up", "Down", "Both", or "All" (default: "Both").

## Value

A named list of character vectors (one per comparison).

## Examples

``` r
v <- example_vista()
#> estimating size factors
#> estimating dispersions
#> gene-wise dispersion estimates
#> mean-dispersion relationship
#> final dispersion estimates
#> fitting model and testing
comp <- names(comparisons(v))[1]
genes <- get_genes_by_regulation(v, sample_comparisons = comp, regulation = 'Up')
str(genes, max.level = 1)
#> List of 1
#>  $ treatment1_VS_control: chr [1:2] "ENSG00000003402" "ENSG00000004799"
```
