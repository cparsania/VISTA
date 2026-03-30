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
if (FALSE) { # \dontrun{
v <- example_vista()
comp <- names(comparisons(v))[1]
if (requireNamespace('org.Mm.eg.db', quietly = TRUE)) {
  out <- try(get_go_enrichment(v, sample_comparison = comp, ont = 'BP', from_type = 'ENSEMBL',
                               orgdb = org.Mm.eg.db::org.Mm.eg.db), silent = TRUE)
  if (!inherits(out, 'try-error')) out
}
} # }
```
