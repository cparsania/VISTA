# Deprecated and defunct arguments in VISTA

VISTA follows the Bioconductor deprecation contract: a renamed argument
warns for one release, becomes defunct in the next, and is removed in
the one after. Deprecated arguments keep working exactly as before while
they warn.

## Value

Not applicable; this page documents argument deprecations only.

## Details

The table below is generated from the package's internal registry, so it
is always in sync with the code. Warnings are classed
`c("vista_deprecated_arg", "deprecatedWarning")` and are emitted at most
once per session per argument, so they will not flood a loop.

To silence them while you migrate:


    suppressWarnings(get_corr_heatmap(v, show_corr_values = TRUE))

To find every deprecated call in your own scripts, promote the warnings
to errors:


    options(warn = 2)

## Inspecting the registry

The authoritative list lives in the package itself rather than being
duplicated here, so it cannot drift out of date. Print it with the
example below, or filter it to the function you are migrating.

## Examples

``` r
# The full table of deprecated arguments and their timelines
VISTA:::.vista_deprecations()[, c("fun", "old_arg", "new_arg", "defunct_in")]
#>                           fun          old_arg          new_arg defunct_in
#> 1            get_corr_heatmap show_corr_values            label      1.4.0
#> 2            get_corr_heatmap  col_corr_values      label_color      1.4.0
#> 3            get_volcano_plot           col_up           colors      1.4.0
#> 4            get_volcano_plot         col_down           colors      1.4.0
#> 5            get_volcano_plot        col_other           colors      1.4.0
#> 6            get_volcano_plot       col_others           colors      1.4.0
#> 7            get_volcano_plot         lab_size       label_size      1.4.0
#> 8                get_pca_plot      sample.seed                       1.4.0
#> 9                get_pca_plot use_vista_colors use_group_colors      1.4.0
#> 10               get_mds_plot use_vista_colors use_group_colors      1.4.0
#> 11              get_umap_plot use_vista_colors use_group_colors      1.4.0
#> 12     get_expression_barplot      facet_scale     facet_scales      1.4.0
#> 13    get_expression_lollipop      facet_scale     facet_scales      1.4.0
#> 14           get_corr_heatmap       cluster_by         order_by      1.4.0
#> 15 get_celltype_group_dotplot            error         errorbar      1.4.0
#> 16      get_deg_count_pieplot            label       label_type      1.4.0
#> 17    get_deg_count_donutplot            label       label_type      1.4.0
#> 18     get_expression_boxplot      comparisons stat_comparisons      1.4.0
#> 19    get_expression_lineplot      comparisons stat_comparisons      1.4.0
#> 20  get_expression_violinplot      comparisons stat_comparisons      1.4.0
#> 21     get_expression_barplot      comparisons stat_comparisons      1.4.0
#> 22               get_pca_plot      top_n_genes            top_n      1.4.0
#> 23               get_mds_plot      top_n_genes            top_n      1.4.0
#> 24              get_umap_plot      top_n_genes            top_n      1.4.0
#> 25    get_expression_lollipop        line_size        linewidth      1.4.0
#> 26    get_foldchange_lollipop        line_size        linewidth      1.4.0
#> 27  get_expression_violinplot  value_transform    log_transform      1.4.0
#> 28    get_expression_lineplot  value_transform    log_transform      1.4.0

# Or just the one function you are migrating
reg <- VISTA:::.vista_deprecations()
reg[reg$fun == "get_corr_heatmap", c("old_arg", "new_arg", "defunct_in")]
#>             old_arg     new_arg defunct_in
#> 1  show_corr_values       label      1.4.0
#> 2   col_corr_values label_color      1.4.0
#> 14       cluster_by    order_by      1.4.0
```
