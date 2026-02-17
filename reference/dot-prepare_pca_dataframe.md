# Prepare PCA result as a tidy data frame

Converts PCA output and metadata into a tidy format suitable for
ggplot2. Includes PC1 and PC2 scores and merges with sample metadata.

## Usage

``` r
.prepare_pca_dataframe(pca, meta)
```

## Arguments

- pca:

  A PCA object returned by
  [`stats::prcomp()`](https://rdrr.io/r/stats/prcomp.html).

- meta:

  A sample metadata data frame, typically from
  [`.prepare_sample_metadata()`](dot-prepare_sample_metadata.md).

## Value

A data frame containing PCA components and sample metadata.
