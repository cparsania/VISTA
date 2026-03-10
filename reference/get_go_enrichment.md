# Run GO enrichment directly from a VISTA comparison

Run GO enrichment directly from a VISTA comparison

## Usage

``` r
get_go_enrichment(
  x,
  sample_comparison,
  regulation = c("Up", "Down", "Both", "All"),
  ont = c("BP", "MF", "CC"),
  from_type = "SYMBOL",
  orgdb = NULL,
  species = "Mus musculus",
  background = NULL,
  ...
)
```

## Arguments

- x:

  A `VISTA` object with DE results.

- sample_comparison:

  Comparison name to use.

- regulation:

  One of `"Up"`, `"Down"`, `"Both"`, or `"All"`; selects genes.

- ont:

  GO ontology: `"BP"`, `"MF"`, or `"CC"`.

- from_type:

  Identifier type in the DE tables (default `"SYMBOL"`).

- orgdb:

  OrgDb object; defaults to mouse/human based on `species`.

- species:

  Species name to infer default OrgDb.

- background:

  Optional background gene set; default uses all features.

- ...:

  Passed to
  [`clusterProfiler::enrichGO()`](https://rdrr.io/pkg/clusterProfiler/man/enrichGO.html).

## Value

A list with `enrich` containing an `enrichResult`.

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
if (requireNamespace('org.Mm.eg.db', quietly = TRUE)) {
  out <- try(get_go_enrichment(v, sample_comparison = comp, ont = 'BP', from_type = 'ENSEMBL',
                               orgdb = org.Mm.eg.db::org.Mm.eg.db), silent = TRUE)
  if (!inherits(out, 'try-error')) out
}
#> --> No gene can be mapped....
#> --> Expected input gene ID: ENSMUSG00000022672,ENSMUSG00000024240,ENSMUSG00000030697,ENSMUSG00000029275,ENSMUSG00000021177,ENSMUSG00000020152
#> --> return NULL...
#> $enrich
#> NULL
#> 
```
