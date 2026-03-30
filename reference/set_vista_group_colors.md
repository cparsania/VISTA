# Set manual group colors in a VISTA object

Updates `metadata(x)$group$colors` using a user-supplied named color
vector. This controls group-level coloring across VISTA plots that use
group mapping (for example PCA, MDS, and expression plots). The color
map must include all groups currently present in the object.

## Usage

``` r
set_vista_group_colors(object, color_map)
```

## Arguments

- object:

  A `VISTA` object.

- color_map:

  Named character vector of colors, with names equal to group labels.

## Value

A modified `VISTA` object with updated group colors.

## Examples

``` r
v <- example_vista()
groups <- unique(as.character(sample_info(v)$cond_long))
gmap <- stats::setNames(c('#1b9e77', '#d95f02')[seq_along(groups)], groups)
set_vista_group_colors(v, gmap)
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
