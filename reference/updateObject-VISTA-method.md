# Update a VISTA object to the current metadata schema

Brings a `VISTA` object created by an older version of the package up to
the current metadata layout. Missing metadata keys are back-filled with
their documented defaults, the schema tag is stamped to the running
version, and the migration is recorded in
`metadata(x)$provenance$updates`.

## Usage

``` r
# S4 method for class 'VISTA'
updateObject(object, ..., verbose = FALSE)
```

## Arguments

- object:

  A `VISTA` object.

- ...:

  Passed to other methods (unused).

- verbose:

  Logical; report what was migrated.

## Value

A `VISTA` object stamped with the current schema version.

## Details

Some content cannot be recovered by migration. In particular, objects
built before VISTA 1.2.0 did not retain raw counts, and normalized
counts are not invertible –
[`counts()`](https://cparsania.github.io/VISTA/reference/counts.md) on
such an object reports that a rebuild is required rather than inventing
values.

Objects carrying a *newer* schema than the running package are left
untouched and reported by
[`validate_vista()`](https://cparsania.github.io/VISTA/reference/validate_vista.md)
as an issue, because newer metadata may carry semantics this version
would misread.

## Examples

``` r
v <- example_vista()

# Simulate an object written by an older release.
S4Vectors::metadata(v)$vista_schema_version <- "0.9.0"
v2 <- BiocGenerics::updateObject(v)
S4Vectors::metadata(v2)$vista_schema_version
#> [1] "1.1.0"
```
