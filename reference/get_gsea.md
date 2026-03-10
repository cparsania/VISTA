# Gene set enrichment analysis (GSEA) from a VISTA comparison

Gene set enrichment analysis (GSEA) from a VISTA comparison

## Usage

``` r
get_gsea(
  x,
  sample_comparison,
  set_type = c("msigdb", "go", "kegg"),
  from_type = "SYMBOL",
  orgdb = NULL,
  species = "Mus musculus",
  msigdb_category = "H",
  msigdb_subcategory = NULL,
  ...
)
```

## Arguments

- x:

  A `VISTA` object with DE results.

- sample_comparison:

  Comparison name to use.

- set_type:

  One of `"msigdb"`, `"go"`, or `"kegg"` selecting the gene set source.

- from_type:

  Identifier type in the DE tables (default `"SYMBOL"`).

- orgdb:

  OrgDb object; defaults to mouse/human based on `species`.

- species:

  Species name to infer default OrgDb.

- msigdb_category, msigdb_subcategory:

  Passed to
  [`msigdbr::msigdbr()`](https://igordot.github.io/msigdbr/reference/msigdbr.html)
  when `set_type = "msigdb"`.

- ...:

  Additional arguments forwarded to the underlying GSEA function:
  [`clusterProfiler::GSEA()`](https://rdrr.io/pkg/clusterProfiler/man/GSEA.html)
  (msigdb TERM2GENE), `gseGO()`, or `gseKEGG()` depending on `set_type`.

## Value

An object returned by this function.

## Examples

``` r
v <- example_vista()
#> estimating size factors
#> estimating dispersions
#> gene-wise dispersion estimates
#> mean-dispersion relationship
#> final dispersion estimates
#> fitting model and testing
comp <- names(comparisons(v))[1]
if (requireNamespace('msigdbr', quietly = TRUE)) {
  out <- try(get_gsea(v, sample_comparison = comp, set_type = 'msigdb', from_type = 'ENSEMBL', species = 'Homo sapiens'), silent = TRUE)
  if (!inherits(out, 'try-error')) out
}
#> 
#> using 'fgsea' for GSEA analysis, please cite Korotkevich et al (2019).
#> preparing geneSet collections...
#> GSEA analysis...
#> no term enriched under specific pvalueCutoff...
#> $enrich
#> #
#> # Gene Set Enrichment Analysis
#> #
#> #...@organism     UNKNOWN 
#> #...@setType      UNKNOWN 
#> #...@geneList     Named num [1:123] 2.918 1.328 1.089 0.984 0.697 ...
#>  - attr(*, "names")= chr [1:123] "ENSG00000004799" "ENSG00000006210" "ENSG00000003402" "ENSG00000003987" ...
#> #...nPerm     
#> #...pvalues adjusted by 'BH' with cutoff <0.05 
#> #...0 enriched terms found
#> 'data.frame':    0 obs. of  8 variables:
#>  $ ID             : chr 
#>  $ Description    : chr 
#>  $ setSize        : int 
#>  $ enrichmentScore: num 
#>  $ NES            : num 
#>  $ pvalue         : num 
#>  $ p.adjust       : num 
#>  $ qvalue         : num 
#> #...Citation
#> S Xu, E Hu, Y Cai, Z Xie, X Luo, L Zhan, W Tang, Q Wang, B Liu, R Wang, W Xie, T Wu, L Xie, G Yu. Using clusterProfiler to characterize multiomics data. Nature Protocols. 2024, 19(11):3292-3320 
#> 
#> 
```
