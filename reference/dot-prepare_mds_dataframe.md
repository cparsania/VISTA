# Prepare MDS Data Frame

Internal helper to convert MDS coordinates and sample metadata into a
plottable data frame.

## Usage

``` r
.prepare_mds_dataframe(mds, meta)
```

## Arguments

- mds:

  A matrix of MDS coordinates from
  [`cmdscale()`](https://rdrr.io/r/stats/cmdscale.html).

- meta:

  A data frame with sample metadata.

## Value

A tidy data frame for MDS plotting.
