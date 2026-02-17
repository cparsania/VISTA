# Run KEGG enrichment directly from a VISTA comparison

Run KEGG enrichment directly from a VISTA comparison

## Usage

``` r
get_kegg_enrichment(
  x,
  sample_comparison,
  regulation = c("Up", "Down", "Both", "All"),
  from_type = "SYMBOL",
  orgdb = NULL,
  species = "Mus musculus",
  kegg_species = NULL,
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

- from_type:

  Identifier type in the DE tables (default `"SYMBOL"`).

- orgdb:

  OrgDb object; defaults to mouse/human based on `species`.

- species:

  Species name to infer default OrgDb.

- kegg_species:

  KEGG organism code (e.g., `"mmu"` or `"hsa"`). If `NULL`, inferred
  from `species`.

- background:

  Optional background gene set; default uses all features.

- ...:

  Passed to
  [`clusterProfiler::enrichGO()`](https://rdrr.io/pkg/clusterProfiler/man/enrichGO.html).
