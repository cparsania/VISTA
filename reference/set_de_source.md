# Set active DE source in a VISTA object

Switches the DE result source used by all downstream VISTA plotting and
accessor functions that read the active metadata slots.

## Usage

``` r
set_de_source(
  object,
  source = c("deseq2", "edger", "limma", "consensus", "active")
)
```

## Arguments

- object:

  A `VISTA` object.

- source:

  One of `"active"`, `"deseq2"`, `"edger"`, `"limma"`, or `"consensus"`.
  When `"active"`, the currently active DE source is kept.

## Value

A modified `VISTA` object with updated active DE source.

## Examples

``` r
v <- example_vista(method = "both")
#> estimating size factors
#> estimating dispersions
#> gene-wise dispersion estimates
#> mean-dispersion relationship
#> final dispersion estimates
#> fitting model and testing
v <- set_de_source(v, "edger")
names(comparisons(v, source = "active"))
#> [1] "treatment1_VS_control"
```
