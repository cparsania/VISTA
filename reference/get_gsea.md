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
NULL
#> NULL
```
