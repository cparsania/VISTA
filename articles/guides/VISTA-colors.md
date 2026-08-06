# Color and Palette Design with VISTA

## Overview

A common challenge in RNA-seq analysis is visual consistency: the same
biological group should keep the same color across PCA, MDS, expression,
and summary plots. When this is not controlled, interpretation becomes
harder and figure legends can become confusing.

This vignette focuses on VISTA’s color system and shows how to:

1.  Discover all built-in palette options.
2.  Keep colors consistent across different plot families.
3.  Change global color style in one place.
4.  Apply manual group/comparison color maps when required.

The floating table of contents on the side provides a quick-scroll panel
for navigating sections.

## Load Packages and Data

``` r

library(VISTA)
library(ggplot2)
library(dplyr)
library(tibble)
library(colorspace)

data("count_data", package = "VISTA")
data("sample_metadata", package = "VISTA")

# Keep vignette runtime moderate while preserving signal
count_small <- count_data[1:1500, ]

# Use all available samples
table(sample_metadata$cond_long)
#> 
#>    control treatment1 
#>          4          4
```

## Build a Baseline VISTA Object

``` r

vista_dark <- create_vista(
  counts = count_small,
  sample_info = sample_metadata,
  column_geneid = "gene_id",
  group_column = "cond_long",
  group_numerator = "treatment1",
  group_denominator = "control",
  group_palette = "Dark 2",
  comparison_palette = "Dark 3",
  min_counts = 5,
  min_replicates = 1
)

vista_dark
#> class: VISTA 
#> dim: 1305 8 
#> metadata(12): de_results de_summary ... design comparison
#> assays(2): norm_counts counts
#> rownames(1305): ENSG00000000003 ENSG00000000419 ... ENSG00000080007
#>   ENSG00000080031
#> rowData names(1): baseMean
#> colnames(8): SRR1039508 SRR1039509 ... SRR1039520 SRR1039521
#> colData names(14): SampleName cell ... sizeFactor sample_names
#> -------- VISTA --------
#> group column: cond_long (control, treatment1)
#> comparisons: treatment1_VS_control
#> DE source: deseq2
#> cutoffs: |log2FC| >= 1, padj <= 0.05
#> raw counts: available via counts()
#> schema: 1.1.0
```

Inspect the stored palette metadata:

``` r

group_palette(vista_dark)
#> [1] "Dark 2"
group_colors(vista_dark)
#>    control treatment1 
#>  "#C87A8A"  "#00A396"
```

## Built-In Palette Choices

VISTA group palettes are based on colorspace qualitative HCL palettes.

``` r

palette_table <- colorspace::hcl_palettes(type = "qualitative") |>
  as.data.frame() |>
  tibble::rownames_to_column("palette") |>
  dplyr::select(palette)

n_palettes <- nrow(palette_table)

cat("Number of built-in qualitative palettes available in VISTA:", n_palettes, "\n\n")
#> Number of built-in qualitative palettes available in VISTA: 9
palette_table
#>    palette
#> 1 Pastel 1
#> 2   Dark 2
#> 3   Dark 3
#> 4    Set 2
#> 5    Set 3
#> 6     Warm
#> 7     Cold
#> 8 Harmonic
#> 9  Dynamic
```

Create a visual palette atlas:

``` r

palette_grid <- do.call(
  rbind,
  lapply(palette_table$palette, function(pal) {
    cols <- colorspace::qualitative_hcl(8, palette = pal)
    data.frame(
      palette = pal,
      idx = seq_along(cols),
      col = cols,
      stringsAsFactors = FALSE
    )
  })
)

palette_grid$palette <- factor(palette_grid$palette, levels = rev(palette_table$palette))

ggplot(palette_grid, aes(x = idx, y = palette, fill = col)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_identity() +
  scale_x_continuous(breaks = 1:8) +
  labs(
    x = "Color index within palette",
    y = NULL,
    title = "VISTA Palette Atlas (Qualitative HCL)"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank())
```

![](VISTA-colors_files/figure-html/palette-atlas-1.png)

## Consistency Across Plot Families

Use the same `VISTA` object and verify that group colors are preserved
across multiple visualization types.

``` r

comp_name <- names(comparisons(vista_dark))[1]
up_genes <- get_genes_by_regulation(
  vista_dark,
  sample_comparisons = comp_name,
  regulation = "Up"
)[[comp_name]]

if (length(up_genes) < 6) {
  genes_demo <- head(rownames(vista_dark), 6)
} else {
  genes_demo <- head(up_genes, 6)
}

genes_demo
#> [1] "ENSG00000003402" "ENSG00000004799" "ENSG00000006788" "ENSG00000008256"
#> [5] "ENSG00000008311" "ENSG00000009413"
```

### PCA (group colors from object metadata)

``` r

get_pca_plot(
  vista_dark,
  label = TRUE,
  shape_by = "cell",
  point_size = 5
)
```

![](VISTA-colors_files/figure-html/pca-consistent-1.png)

### MDS (same group colors)

``` r

get_mds_plot(
  vista_dark,
  label = TRUE,
  shape_by = "cell",
  point_size = 5
)
```

![](VISTA-colors_files/figure-html/mds-consistent-1.png)

### Expression barplot (same group colors)

``` r

get_expression_barplot(
  vista_dark,
  genes = genes_demo[1:4],
  log_transform = TRUE
)
```

![](VISTA-colors_files/figure-html/expr-bar-consistent-1.png)

### Expression boxplot (same group colors)

``` r

get_expression_boxplot(
  vista_dark,
  genes = genes_demo[1:4],
  log_transform = TRUE,
  pool_genes = FALSE,
  by = "gene",
  facet_by = "none",
  fill_by = "group"
)
```

![](VISTA-colors_files/figure-html/expr-box-consistent-1.png)

## Switching Global Style in One Argument

Rebuild with a different global palette. No plot-specific color edits
are needed.

``` r

vista_warm <- create_vista(
  counts = count_small,
  sample_info = sample_metadata,
  column_geneid = "gene_id",
  group_column = "cond_long",
  group_numerator = "treatment1",
  group_denominator = "control",
  group_palette = "Warm",
  comparison_palette = "Set 3",
  min_counts = 5,
  min_replicates = 1
)

group_palette(vista_warm)
#> [1] "Warm"
group_colors(vista_warm)
#>    control treatment1 
#>  "#ABB065"  "#E093C3"
```

``` r

p_dark <- get_pca_plot(vista_dark, label = FALSE, point_size = 5) +
  ggtitle("Dark 2")

p_warm <- get_pca_plot(vista_warm, label = FALSE, point_size = 5) +
  ggtitle("Warm")

if (requireNamespace("patchwork", quietly = TRUE)) {
  patchwork::wrap_plots(p_dark, p_warm, ncol = 2)
} else {
  p_dark
  p_warm
}
```

![](VISTA-colors_files/figure-html/pca-palette-compare-1.png)

## Comparison-Level Colors

For multi-comparison workflows, VISTA also stores comparison colors.

``` r

# Build two directed comparisons to demonstrate comparison color mapping
vista_multi <- create_vista(
  counts = count_small,
  sample_info = sample_metadata,
  column_geneid = "gene_id",
  group_column = "cond_long",
  group_numerator = c("treatment1", "control"),
  group_denominator = c("control", "treatment1"),
  group_palette = "Dark 2",
  comparison_palette = "Set 3",
  min_counts = 5,
  min_replicates = 1
)

S4Vectors::metadata(vista_multi)$comparison$colors
#> treatment1_VS_control control_VS_treatment1 
#>             "#FFB3B5"             "#61D8D6"
```

``` r

get_foldchange_barplot(
  vista_multi,
  genes = head(rownames(vista_multi), 6),
  facet_by = "none"
)
```

![](VISTA-colors_files/figure-html/fc-bar-comparison-colors-1.png)

## Manual Color Maps (Advanced)

If your lab or journal has fixed brand colors, you can set explicit
named vectors using built-in VISTA setters.

Available manual setters:

- `set_vista_group_colors(vista, color_map)`
- `set_vista_comparison_colors(vista, color_map)`

``` r

custom_group_cols <- c(
  control = "#264653",
  treatment1 = "#E76F51"
)

custom_comp_cols <- c(
  treatment1_VS_control = "#6C5CE7",
  control_VS_treatment1 = "#00A896"
)

vista_custom <- set_vista_group_colors(vista_dark, custom_group_cols)
vista_multi_custom <- set_vista_comparison_colors(vista_multi, custom_comp_cols)

group_colors(vista_custom)
#>    control treatment1 
#>  "#264653"  "#E76F51"
S4Vectors::metadata(vista_multi_custom)$comparison$colors
#> treatment1_VS_control control_VS_treatment1 
#>             "#6C5CE7"             "#00A896"
```

``` r

p_custom_group <- get_pca_plot(
  vista_custom,
  label = FALSE,
  point_size = 5
) + ggtitle("Custom group colors")

p_custom_comp <- get_foldchange_barplot(
  vista_multi_custom,
  genes = head(rownames(vista_multi_custom), 6),
  facet_by = "none"
) + ggtitle("Custom comparison colors")

if (requireNamespace("patchwork", quietly = TRUE)) {
  patchwork::wrap_plots(p_custom_group, p_custom_comp, ncol = 2)
} else {
  p_custom_group
  p_custom_comp
}
```

![](VISTA-colors_files/figure-html/custom-color-plots-1.png)

## Practical Guidance

- Use `group_palette` and `comparison_palette` in
  [`create_vista()`](https://cparsania.github.io/VISTA/reference/create_vista.md)
  for fast, global styling.
- Use `group_colors(vista)` to audit mapping before generating
  publication panels.
- Keep one palette per manuscript/project to preserve interpretation
  consistency.
- Switch to manual named vectors only when you need strict journal or
  brand colors.

## Session Information

``` r

sessionInfo()
#> R version 4.6.1 (2026-06-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] colorspace_2.1-3 tibble_3.3.1     dplyr_1.2.1      ggplot2_4.0.3   
#> [5] VISTA_1.1.5      BiocStyle_2.40.0
#> 
#> loaded via a namespace (and not attached):
#>   [1] RColorBrewer_1.1-3          jsonlite_2.0.0             
#>   [3] tidydr_0.0.6                magrittr_2.0.5             
#>   [5] ggtangle_0.1.2              farver_2.1.2               
#>   [7] rmarkdown_2.31              fs_2.1.0                   
#>   [9] ragg_1.5.2                  vctrs_0.7.3                
#>  [11] memoise_2.0.1               ggtree_4.2.0               
#>  [13] rstatix_1.1.0               htmltools_0.5.9            
#>  [15] S4Arrays_1.12.0             curl_7.1.0                 
#>  [17] broom_1.0.13                Formula_1.2-6              
#>  [19] SparseArray_1.12.2          gridGraphics_0.5-1         
#>  [21] sass_0.4.10                 bslib_0.12.0               
#>  [23] htmlwidgets_1.6.4           desc_1.4.3                 
#>  [25] plyr_1.8.9                  httr2_1.3.0                
#>  [27] cachem_1.1.0                igraph_2.3.3               
#>  [29] lifecycle_1.0.5             pkgconfig_2.0.3            
#>  [31] gson_0.2.1                  Matrix_1.7-5               
#>  [33] R6_2.6.1                    fastmap_1.2.0              
#>  [35] MatrixGenerics_1.24.0       digest_0.6.39              
#>  [37] aplot_0.3.1                 enrichplot_1.32.0          
#>  [39] ggnewscale_0.5.2            GGally_2.4.0               
#>  [41] patchwork_1.3.2             AnnotationDbi_1.74.0       
#>  [43] S4Vectors_0.50.1            aisdk_1.4.12               
#>  [45] ps_1.9.3                    DESeq2_1.52.0              
#>  [47] textshaping_1.0.5           GenomicRanges_1.64.0       
#>  [49] RSQLite_3.53.3              ggpubr_1.0.0               
#>  [51] labeling_0.4.3              polyclip_1.10-7            
#>  [53] httr_1.4.8                  abind_1.4-8                
#>  [55] compiler_4.6.1              withr_3.0.3                
#>  [57] bit64_4.8.2                 fontquiver_0.2.1           
#>  [59] backports_1.5.1             S7_0.2.2                   
#>  [61] BiocParallel_1.46.0         carData_3.0-6              
#>  [63] DBI_1.3.0                   ggstats_0.13.0             
#>  [65] ggforce_0.5.0               ggsignif_0.6.4             
#>  [67] MASS_7.3-65                 rappdirs_0.3.4             
#>  [69] DelayedArray_0.38.2         tools_4.6.1                
#>  [71] otel_0.2.0                  scatterpie_0.2.6           
#>  [73] ape_5.8-1                   msigdbr_26.1.0             
#>  [75] glue_1.8.1                  callr_3.8.0                
#>  [77] nlme_3.1-169                GOSemSim_2.38.3            
#>  [79] grid_4.6.1                  cluster_2.1.8.2            
#>  [81] reshape2_1.4.5              generics_0.1.4             
#>  [83] gtable_0.3.6                tidyr_1.3.2                
#>  [85] car_3.1-5                   XVector_0.52.0             
#>  [87] BiocGenerics_0.58.1         ggrepel_0.9.8              
#>  [89] pillar_1.11.1               stringr_1.6.0              
#>  [91] babelgene_22.9              limma_3.68.4               
#>  [93] yulab.utils_0.2.4           splines_4.6.1              
#>  [95] tweenr_2.0.3                treeio_1.36.1              
#>  [97] lattice_0.22-9              bit_4.6.0                  
#>  [99] tidyselect_1.2.1            fontLiberation_0.1.0       
#> [101] GO.db_3.23.1                locfit_1.5-9.12            
#> [103] Biostrings_2.80.1           knitr_1.51                 
#> [105] fontBitstreamVera_0.1.1     bookdown_0.47              
#> [107] IRanges_2.46.0              Seqinfo_1.2.0              
#> [109] edgeR_4.10.1                SummarizedExperiment_1.42.0
#> [111] stats4_4.6.1                xfun_0.60                  
#> [113] Biobase_2.72.0              statmod_1.5.2              
#> [115] matrixStats_1.5.0           stringi_1.8.9              
#> [117] lazyeval_0.2.3              ggfun_0.2.1                
#> [119] yaml_2.3.12                 evaluate_1.0.5             
#> [121] codetools_0.2-20            qvalue_2.44.0              
#> [123] gdtools_0.5.1               BiocManager_1.30.27        
#> [125] ggplotify_0.1.3             cli_3.6.6                  
#> [127] systemfonts_1.3.2           processx_3.9.0             
#> [129] jquerylib_0.1.4             Rcpp_1.1.2                 
#> [131] png_0.1-9                   parallel_4.6.1             
#> [133] assertthat_0.2.1            pkgdown_2.2.1              
#> [135] blob_1.3.0                  clusterProfiler_4.20.0     
#> [137] DOSE_4.6.0                  tidytree_0.4.8             
#> [139] ggiraph_0.9.6               enrichit_0.2.1             
#> [141] scales_1.4.0                purrr_1.2.2                
#> [143] crayon_1.5.3                rlang_1.3.0                
#> [145] KEGGREST_1.52.2
```
