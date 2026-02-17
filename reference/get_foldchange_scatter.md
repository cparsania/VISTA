# Fold-change scatterplot between two comparisons

Plots log2 fold changes from two stored comparisons against each other,
with points coloured by concordant/discordant regulation based on the
cutoffs saved in the `VISTA` object.

## Usage

``` r
get_foldchange_scatter(
  x,
  sample_comparisons,
  label_top_n = 0,
  point_alpha = 0.5,
  use_hex = FALSE,
  method = c("pearson", "spearman")
)
```

## Arguments

- x:

  A `VISTA` object containing DE results.

- sample_comparisons:

  Character vector of length 2 naming the comparisons.

- label_top_n:

  Integer; number of most extreme points to label (by \|log2FC1\| +
  \|log2FC2\|).

- point_alpha:

  Point transparency.

- use_hex:

  Logical; use hex binning instead of points for dense datasets.

- method:

  Correlation method for the subtitle; `"pearson"` or `"spearman"`.

## Value

A `ggplot2` object.

## Details

Points are coloured by concordance status using fixed colours:

- Up/Up = `#1b9e77`

- Down/Down = `#7570b3`

- Up/Down = `#d95f02`

- Down/Up = `#e7298a`

- Other = `grey70`

Regulation is derived from the `log2fc` and `pval` cutoffs stored in
`cutoffs(x)` (and `p_value_type` from the same list, defaulting to
`"padj"`).
