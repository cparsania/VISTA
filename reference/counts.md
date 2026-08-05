# Raw counts from a VISTA object

Returns the filtered raw count matrix stored alongside the normalized
counts. Objects built before VISTA 1.2.0, or with
`create_vista(keep_raw_counts = FALSE)`, do not carry one.

## Usage

``` r
# S4 method for class 'VISTA'
counts(object, ...)
```

## Arguments

- object:

  A `VISTA` object.

- ...:

  Ignored.

## Value

A numeric matrix of raw counts with genes in rows and samples in
columns.

## Details

Normalized counts are not invertible, so
[`BiocGenerics::updateObject()`](https://rdrr.io/pkg/BiocGenerics/man/updateObject.html)
cannot back-fill this assay for an older object; it has to be rebuilt
with
[`create_vista()`](https://cparsania.github.io/VISTA/reference/create_vista.md).

## Examples

``` r
v <- example_vista()
counts(v)[seq_len(3), seq_len(3)]
#>                 SRR1039508 SRR1039512 SRR1039516
#> ENSG00000000003        679        873       1138
#> ENSG00000000419        467        621        587
#> ENSG00000000457        260        263        245
```
