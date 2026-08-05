# Preparing Counts and Metadata for VISTA

## Overview

VISTA expects two main inputs for
[`create_vista()`](https://cparsania.github.io/VISTA/reference/create_vista.md):

1.  A count table with one gene identifier column plus one column per
    sample.
2.  A sample metadata table (`sample_info`) containing `sample_names`
    and the grouping variables used in the analysis.

In practice, users often start from count tables whose column names are
file paths, alignment outputs, or sequencing-derived sample identifiers.
Metadata can also be incomplete or missing entirely. VISTA’s input
preparation helpers are designed to standardize these objects before
differential analysis starts:

- [`read_vista_counts()`](https://cparsania.github.io/VISTA/reference/read_vista_counts.md)
- [`derive_vista_metadata()`](https://cparsania.github.io/VISTA/reference/derive_vista_metadata.md)
- [`read_vista_metadata()`](https://cparsania.github.io/VISTA/reference/read_vista_metadata.md)
- [`match_vista_inputs()`](https://cparsania.github.io/VISTA/reference/match_vista_inputs.md)

## Load Example Data

``` r

library(VISTA)
library(dplyr)
library(tibble)

data("count_data", package = "VISTA")
data("sample_metadata", package = "VISTA")

dim(count_data)
#> [1] 63677     9
head(sample_metadata[, c("sample_names", "cond_long")])
#> # A tibble: 6 × 2
#>   sample_names cond_long 
#>   <chr>        <chr>     
#> 1 SRR1039508   control   
#> 2 SRR1039509   treatment1
#> 3 SRR1039512   control   
#> 4 SRR1039513   treatment1
#> 5 SRR1039516   control   
#> 6 SRR1039517   treatment1
```

## 1. Standardize Counts

For an ordinary in-memory count table,
[`read_vista_counts()`](https://cparsania.github.io/VISTA/reference/read_vista_counts.md)
keeps the existing API contract and returns a
[`create_vista()`](https://cparsania.github.io/VISTA/reference/create_vista.md)-ready
structure.

``` r

prepared_counts <- read_vista_counts(
  count_data[seq_len(50), ],
  format = "matrix",
  gene_id_column = "gene_id"
)

names(prepared_counts)
#> [1] "counts"          "row_data"        "column_geneid"   "sample_names"   
#> [5] "sample_name_map" "input_format"    "report"
head(prepared_counts$counts[, seq_len(4)])
#>           gene_id SRR1039508 SRR1039509 SRR1039512
#> 1 ENSG00000000003        679        448        873
#> 2 ENSG00000000005          0          0          0
#> 3 ENSG00000000419        467        515        621
#> 4 ENSG00000000457        260        211        263
#> 5 ENSG00000000460         60         55         40
#> 6 ENSG00000000938          0          0          2
head(prepared_counts$row_data)
#>                         gene_id
#> ENSG00000000003 ENSG00000000003
#> ENSG00000000005 ENSG00000000005
#> ENSG00000000419 ENSG00000000419
#> ENSG00000000457 ENSG00000000457
#> ENSG00000000460 ENSG00000000460
#> ENSG00000000938 ENSG00000000938
```

The returned object includes:

- `counts`: standardized count table with `gene_id`
- `row_data`: feature metadata aligned to the counts
- `sample_names`: standardized sample columns
- `sample_name_map`: mapping from original to repaired sample names

## 2. Repair File-Derived Sample Names

RNA-seq count tables often inherit sample columns from alignment or
quantification files. VISTA can repair these names conservatively when
they are clearly file-derived and the repaired names remain unique.

``` r

counts_paths <- count_data[seq_len(20), c("gene_id", sample_metadata$sample_names[seq_len(4)]), drop = FALSE]

colnames(counts_paths)[2:5] <- c(
  "/proj/run/03_alignment/HET_1_U/HET_1_U_star_alignAligned.sortedByCoord.out.bam",
  "/proj/run/03_alignment/HET_1_ovary/HET_1_ovary_star_alignAligned.sortedByCoord.out.bam",
  "/proj/run/03_alignment/WT_1_U/WT_1_U_star_alignAligned.sortedByCoord.out.bam",
  "/proj/run/03_alignment/WT_1_ovary/WT_1_ovary_star_alignAligned.sortedByCoord.out.bam"
)

prepared_paths <- read_vista_counts(
  counts_paths,
  format = "matrix",
  gene_id_column = "gene_id",
  repair_sample_names = "auto"
)

prepared_paths$sample_name_map
#>                                                                                 original
#> 1         /proj/run/03_alignment/HET_1_U/HET_1_U_star_alignAligned.sortedByCoord.out.bam
#> 2 /proj/run/03_alignment/HET_1_ovary/HET_1_ovary_star_alignAligned.sortedByCoord.out.bam
#> 3           /proj/run/03_alignment/WT_1_U/WT_1_U_star_alignAligned.sortedByCoord.out.bam
#> 4   /proj/run/03_alignment/WT_1_ovary/WT_1_ovary_star_alignAligned.sortedByCoord.out.bam
#>      repaired
#> 1     HET_1_U
#> 2 HET_1_ovary
#> 3      WT_1_U
#> 4  WT_1_ovary
```

Automatic repair currently handles patterns such as:

- full file paths
- STAR outputs (`Aligned.sortedByCoord.out.bam`, `ReadsPerGene.out.tab`)
- RSEM outputs (`.genes.results`, `.isoforms.results`)
- generic quantification files such as `quant.sf`
- lane/read suffixes like `_S1_L001_R1_001`

## 3. Derive Starter Metadata from Sample Names

If you do not yet have a metadata sheet,
[`derive_vista_metadata()`](https://cparsania.github.io/VISTA/reference/derive_vista_metadata.md)
can create a starter `sample_info` table from the count sample names.

### Split-based parsing

``` r

derived_split <- derive_vista_metadata(
  prepared_paths,
  parser = "split",
  split = "_",
  fields = c("genotype", "replicate", "tissue")
)

derived_split
#>   sample_names genotype replicate tissue
#> 1      HET_1_U      HET         1      U
#> 2  HET_1_ovary      HET         1  ovary
#> 3       WT_1_U       WT         1      U
#> 4   WT_1_ovary       WT         1  ovary
```

### Regex-based parsing

For sequencing accessions or other structured labels, regex parsing
gives more control:

``` r

derived_regex <- derive_vista_metadata(
  counts = NULL,
  sample_names = c("SRR1039508", "SRR1039509", "SRR1039512"),
  parser = "regex",
  pattern = "SRR(\\d+)",
  fields = "run_id"
)

derived_regex
#>   sample_names  run_id
#> 1   SRR1039508 1039508
#> 2   SRR1039509 1039509
#> 3   SRR1039512 1039512
```

### Template mode

Template mode adds empty placeholders for columns users commonly fill in
next, such as `group` and `batch`.

``` r

derive_vista_metadata(
  counts = NULL,
  sample_names = c("sampleA", "sampleB"),
  parser = "none",
  return_type = "template"
)
#>   sample_names group batch
#> 1      sampleA  <NA>  <NA>
#> 2      sampleB  <NA>  <NA>
```

## 4. Standardize Existing Metadata

If a metadata table already exists,
[`read_vista_metadata()`](https://cparsania.github.io/VISTA/reference/read_vista_metadata.md)
standardizes it into the form expected by
[`create_vista()`](https://cparsania.github.io/VISTA/reference/create_vista.md).

``` r

prepared_samples <- read_vista_metadata(sample_metadata)

head(prepared_samples[, c("sample_names", "cond_long")])
#>            sample_names  cond_long
#> SRR1039508   SRR1039508    control
#> SRR1039509   SRR1039509 treatment1
#> SRR1039512   SRR1039512    control
#> SRR1039513   SRR1039513 treatment1
#> SRR1039516   SRR1039516    control
#> SRR1039517   SRR1039517 treatment1
```

VISTA will infer `sample_names` from:

- an existing `sample_names` column
- non-default rownames
- common aliases such as `sample`, `sample_id`, or `Run`

## 5. Match Counts and Metadata

[`match_vista_inputs()`](https://cparsania.github.io/VISTA/reference/match_vista_inputs.md)
aligns count columns and metadata rows so they can be passed directly
into
[`create_vista()`](https://cparsania.github.io/VISTA/reference/create_vista.md).

``` r

matched_inputs <- match_vista_inputs(prepared_counts, prepared_samples)

matched_inputs$column_geneid
#> [1] "gene_id"
identical(matched_inputs$sample_info$sample_names, colnames(matched_inputs$counts)[-1])
#> [1] TRUE
```

When appropriate, you can drop unmatched samples instead of erroring:

``` r

sample_extra <- prepared_samples[seq_len(4), , drop = FALSE]
extra_row <- sample_extra[1, , drop = FALSE]
extra_row$sample_names <- "extra_sample"
if ("cond_long" %in% colnames(extra_row)) extra_row$cond_long <- "control"
if ("dex" %in% colnames(extra_row)) extra_row$dex <- "untrt"
if ("cell" %in% colnames(extra_row)) extra_row$cell <- "N999"
rownames(extra_row) <- "extra_sample"
sample_extra <- rbind(sample_extra, extra_row)

matched_drop <- match_vista_inputs(
  prepared_counts,
  sample_extra,
  drop_unmatched = TRUE
)

matched_drop$report
#> $n_genes
#> [1] 50
#> 
#> $n_samples
#> [1] 4
#> 
#> $dropped_from_counts
#> [1] "SRR1039516" "SRR1039517" "SRR1039520" "SRR1039521"
#> 
#> $dropped_from_sample_info
#> [1] "extra_sample"
```

## 6. Create a VISTA Object

Once inputs are standardized,
[`create_vista()`](https://cparsania.github.io/VISTA/reference/create_vista.md)
is unchanged.

``` r

vista <- create_vista(
  counts = matched_inputs$counts,
  sample_info = matched_inputs$sample_info,
  column_geneid = matched_inputs$column_geneid,
  group_column = "cond_long",
  group_numerator = "treatment1",
  group_denominator = "control",
  method = "deseq2",
  min_counts = 5,
  min_replicates = 1
)

vista
#> class: VISTA 
#> dim: 44 8 
#> metadata(12): de_results de_summary ... design comparison
#> assays(2): norm_counts counts
#> rownames(44): ENSG00000000003 ENSG00000000419 ... ENSG00000004399
#>   ENSG00000004455
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

## Practical Guidance

- Use
  [`read_vista_counts()`](https://cparsania.github.io/VISTA/reference/read_vista_counts.md)
  first whenever count columns come from external tools or file names.
- Use
  [`derive_vista_metadata()`](https://cparsania.github.io/VISTA/reference/derive_vista_metadata.md)
  when you only have sample names and need a reviewable starter metadata
  sheet.
- Use
  [`read_vista_metadata()`](https://cparsania.github.io/VISTA/reference/read_vista_metadata.md)
  when you already have a metadata table and want VISTA to standardize
  it.
- Use
  [`match_vista_inputs()`](https://cparsania.github.io/VISTA/reference/match_vista_inputs.md)
  before
  [`create_vista()`](https://cparsania.github.io/VISTA/reference/create_vista.md)
  to catch alignment issues early.

## Session Info

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
#> [1] tibble_3.3.1     dplyr_1.2.1      VISTA_1.1.3      BiocStyle_2.40.0
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
#>  [39] colorspace_2.1-3            ggnewscale_0.5.2           
#>  [41] GGally_2.4.0                patchwork_1.3.2            
#>  [43] AnnotationDbi_1.74.0        S4Vectors_0.50.1           
#>  [45] aisdk_1.4.12                ps_1.9.3                   
#>  [47] DESeq2_1.52.0               textshaping_1.0.5          
#>  [49] GenomicRanges_1.64.0        RSQLite_3.53.3             
#>  [51] ggpubr_1.0.0                polyclip_1.10-7            
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
#>  [85] utf8_1.2.6                  car_3.1-5                  
#>  [87] XVector_0.52.0              BiocGenerics_0.58.1        
#>  [89] ggrepel_0.9.8               pillar_1.11.1              
#>  [91] stringr_1.6.0               babelgene_22.9             
#>  [93] limma_3.68.4                yulab.utils_0.2.4          
#>  [95] splines_4.6.1               tweenr_2.0.3               
#>  [97] treeio_1.36.1               lattice_0.22-9             
#>  [99] bit_4.6.0                   tidyselect_1.2.1           
#> [101] fontLiberation_0.1.0        GO.db_3.23.1               
#> [103] locfit_1.5-9.12             Biostrings_2.80.1          
#> [105] knitr_1.51                  fontBitstreamVera_0.1.1    
#> [107] bookdown_0.47               IRanges_2.46.0             
#> [109] Seqinfo_1.2.0               edgeR_4.10.1               
#> [111] SummarizedExperiment_1.42.0 stats4_4.6.1               
#> [113] xfun_0.60                   Biobase_2.72.0             
#> [115] statmod_1.5.2               matrixStats_1.5.0          
#> [117] stringi_1.8.9               lazyeval_0.2.3             
#> [119] ggfun_0.2.1                 yaml_2.3.12                
#> [121] evaluate_1.0.5              codetools_0.2-20           
#> [123] qvalue_2.44.0               gdtools_0.5.1              
#> [125] BiocManager_1.30.27         ggplotify_0.1.3            
#> [127] cli_3.6.6                   systemfonts_1.3.2          
#> [129] processx_3.9.0              jquerylib_0.1.4            
#> [131] Rcpp_1.1.2                  png_0.1-9                  
#> [133] parallel_4.6.1              assertthat_0.2.1           
#> [135] pkgdown_2.2.1               ggplot2_4.0.3              
#> [137] blob_1.3.0                  clusterProfiler_4.20.0     
#> [139] DOSE_4.6.0                  tidytree_0.4.8             
#> [141] ggiraph_0.9.6               enrichit_0.2.1             
#> [143] scales_1.4.0                purrr_1.2.2                
#> [145] crayon_1.5.3                rlang_1.3.0                
#> [147] KEGGREST_1.52.2
```
