# Set active DE source in a VISTA object

Switches the DE result source used by all downstream VISTA plotting and
accessor functions that read the active metadata slots.

## Usage

``` r
set_de_source(object, source = c("deseq2", "edger", "limma", "consensus"))
```

## Arguments

- object:

  A `VISTA` object.

- source:

  One of `"deseq2"`, `"edger"`, `"limma"`, or `"consensus"`.

## Value

A modified `VISTA` object with updated active DE source.

## Examples

``` r
NULL
#> NULL
```
