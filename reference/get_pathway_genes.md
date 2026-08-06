# Extract genes from enriched pathways

Parses pathway-level gene members from an enrichment result and returns
either a long table, pathway-indexed list, or a unique vector of genes.

## Usage

``` r
get_pathway_genes(
  x,
  pathways = NULL,
  top_n = 10,
  pathway_column = c("Description", "ID"),
  gene_column = c("auto", "geneID", "core_enrichment"),
  gene_sep = "/",
  return_type = c("long", "list", "vector")
)
```

## Arguments

- x:

  An `enrichResult`/`gseaResult`, or a list containing element `enrich`
  (e.g. output from
  [`get_msigdb_enrichment()`](https://cparsania.github.io/VISTA/reference/get_msigdb_enrichment.md)).

- pathways:

  Optional character vector of pathway names to keep. Matches against
  `pathway_column`.

- top_n:

  Number of top pathways to use when `pathways` is `NULL`. Ranking uses
  `p.adjust` (then `pvalue`) when available. Default: `10`.

- pathway_column:

  Which enrichment column to match pathway names against:
  `"Description"` (default) or `"ID"`.

- gene_column:

  Which column stores pathway members. `"auto"` (default) uses
  `"geneID"` when present, otherwise `"core_enrichment"`.

- gene_sep:

  Delimiter used in pathway gene strings (default `"/"`).

- return_type:

  One of `"long"`, `"list"`, or `"vector"`.

## Value

Depending on `return_type`:

- `"long"`: data frame with `pathway_id`, `pathway`, and `gene`.

- `"list"`: named list of character vectors (genes per pathway).

- `"vector"`: unique character vector of genes.

## Examples

``` r
if (requireNamespace("msigdbr", quietly = TRUE)) {
  vista <- example_vista()
  msig <- get_msigdb_enrichment(
    vista,
    sample_comparison = names(comparisons(vista))[1],
    regulation = "Both",
    msigdb_category = "H",
    from_type = "ENSEMBL"
  )
  pathway_tbl <- get_pathway_genes(msig, top_n = 5, return_type = "long")
  head(pathway_tbl)
}
#> Warning: qvalue::qvalue() failed, returning NA for qvalue. Error: missing values and NaN's not allowed if 'na.rm' is FALSE
#>     pathway_id      pathway            gene
#> 1 UN_ANNOTATED UN_ANNOTATED ENSG00000003402
#> 2 UN_ANNOTATED UN_ANNOTATED ENSG00000004799

# \donttest{
data("count_data", package = "VISTA")
data("sample_metadata", package = "VISTA")

vista <- create_vista(
  counts = count_data,
  sample_info = sample_metadata,
  column_geneid = "gene_id",
  group_column = "cond_long",
  group_numerator = "treatment1",
  group_denominator = "control"
)
#> estimating size factors
#> estimating dispersions
#> gene-wise dispersion estimates
#> mean-dispersion relationship
#> final dispersion estimates
#> fitting model and testing

msig <- get_msigdb_enrichment(
  vista,
  sample_comparison = names(comparisons(vista))[1],
  regulation = "Up",
  species = "Homo sapiens",
  from_type = "ENSEMBL"
)
#> 

pathway_tbl <- get_pathway_genes(msig, top_n = 5, return_type = "long")
head(pathway_tbl)
#>                         pathway_id                          pathway
#> 1 HALLMARK_TNFA_SIGNALING_VIA_NFKB HALLMARK_TNFA_SIGNALING_VIA_NFKB
#> 2 HALLMARK_TNFA_SIGNALING_VIA_NFKB HALLMARK_TNFA_SIGNALING_VIA_NFKB
#> 3 HALLMARK_TNFA_SIGNALING_VIA_NFKB HALLMARK_TNFA_SIGNALING_VIA_NFKB
#> 4 HALLMARK_TNFA_SIGNALING_VIA_NFKB HALLMARK_TNFA_SIGNALING_VIA_NFKB
#> 5 HALLMARK_TNFA_SIGNALING_VIA_NFKB HALLMARK_TNFA_SIGNALING_VIA_NFKB
#> 6 HALLMARK_TNFA_SIGNALING_VIA_NFKB HALLMARK_TNFA_SIGNALING_VIA_NFKB
#>              gene
#> 1 ENSG00000143878
#> 2 ENSG00000122641
#> 3 ENSG00000067082
#> 4 ENSG00000162616
#> 5 ENSG00000179094
#> 6 ENSG00000159200
# }
```
