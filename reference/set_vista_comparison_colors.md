# Set manual comparison colors in a VISTA object

Updates `metadata(x)$comparison$colors` using a user-supplied named
color vector. Comparison colors are used in fold-change and
comparison-level plots. The color map must include all currently active
comparisons.

## Usage

``` r
set_vista_comparison_colors(object, color_map)
```

## Arguments

- object:

  A `VISTA` object.

- color_map:

  Named character vector of colors, with names equal to comparison names
  (for example `"A_VS_B"`).

## Value

A modified `VISTA` object with updated comparison colors.

## Examples

``` r
v <- example_vista()
comps <- names(comparisons(v))
if (length(comps)) {
  cmap <- stats::setNames(rep('#1b9e77', length(comps)), comps)
  set_vista_comparison_colors(v, cmap)
}
#> class: SummarizedExperiment 
#> dim: 123 6 
#> metadata(12): de_results de_summary ... design comparison
#> assays(1): norm_counts
#> rownames(123): ENSG00000000003 ENSG00000000419 ... ENSG00000006607
#>   ENSG00000006625
#> rowData names(1): baseMean
#> colnames(6): SRR1039508 SRR1039512 ... SRR1039513 SRR1039517
#> colData names(14): SampleName cell ... sizeFactor sample_names
```
