# Prepare sample metadata with optional filtering and group ordering

This internal helper returns a sample metadata `data.frame` from a VISTA
object, optionally filtered by specific sample groups. The function also
ensures the grouping column is treated as a factor with appropriate
level ordering, either based on user input (`sample_group`) or on
original appearance.

## Usage

``` r
.prepare_sample_metadata(x, sample_group = NULL, group_column = NULL)
```

## Arguments

- x:

  A VISTA object.

- sample_group:

  Optional character vector of groups to include, based on the
  `group_column`.

## Value

A filtered `data.frame` containing sample metadata and a `sample`
column.

## Details

This is a general-purpose metadata preparation function intended to
support multiple downstream plotting or reporting utilities.
