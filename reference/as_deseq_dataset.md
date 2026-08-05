# Convert a VISTA object to a DESeqDataSet

Builds a
[`DESeq2::DESeqDataSet`](https://rdrr.io/pkg/DESeq2/man/DESeqDataSet.html)
from the raw counts retained in a `VISTA` object, so an analysis can be
continued or re-run with DESeq2 directly.

## Usage

``` r
as_deseq_dataset(x, design = NULL)
```

## Arguments

- x:

  A `VISTA` object carrying a raw `counts` assay.

- design:

  A model formula. Defaults to `~ <group column>` using the grouping
  column stored in the object.

## Value

A `DESeqDataSet`.

## Details

This is deliberately explicit rather than making
`DESeq2::DESeqDataSet(v, design)` work implicitly: the inherited
`SummarizedExperiment` method would pick up whichever assay comes first,
which is `norm_counts`, and fail because those values are not integers.

Counts are rounded to integers, as DESeq2 requires.

## Examples

``` r
v <- example_vista()
if (requireNamespace("DESeq2", quietly = TRUE)) {
  dds <- as_deseq_dataset(v)
  class(dds)
}
#> converting counts to integer mode
#> [1] "DESeqDataSet"
#> attr(,"package")
#> [1] "DESeq2"
```
