# Chord diagram of enrichment gene–pathway relationships

Draws a chord diagram linking genes to the enriched pathways they belong
to. Chords can be coloured by fold-change, regulation direction, or
source pathway.

## Usage

``` r
get_enrichment_chord(
  x,
  vista_obj = NULL,
  comparison = NULL,
  pathways = NULL,
  top_n = 8,
  pathway_column = c("Description", "ID"),
  gene_column = c("auto", "geneID", "core_enrichment"),
  gene_sep = "/",
  min_pathways = 1,
  max_genes = 50,
  gene_order_by = c("none", "foldchange", "abs_foldchange"),
  gene_id_column = NULL,
  display_id = NULL,
  color_by = c("foldchange", "regulation", "pathway"),
  up_color = "#D73027",
  down_color = "#1A9850",
  other_color = "grey70",
  pathway_colors = NULL,
  transparency = 0.4,
  gap_degree = 2,
  label_cex = 0.7,
  title = NULL
)
```

## Arguments

- x:

  An `enrichResult`, `gseaResult`, or `compareClusterResult` from
  clusterProfiler, or a list containing an `enrich` element (e.g. output
  of [`get_msigdb_enrichment()`](get_msigdb_enrichment.md)).

- vista_obj:

  Optional `VISTA` object. Required when `color_by` is `"foldchange"` or
  `"regulation"`.

- comparison:

  Character scalar naming the DE comparison in `vista_obj` to pull
  log2FC values from. Required when `vista_obj` is supplied.

- pathways:

  Optional character vector of pathway names to include. Matches against
  `pathway_column`.

- top_n:

  Number of top pathways to display when `pathways` is `NULL` (default
  `8`).

- pathway_column:

  Column in the enrichment result to match pathway names:
  `"Description"` (default) or `"ID"`.

- gene_column:

  Column storing gene members: `"auto"` (default), `"geneID"`, or
  `"core_enrichment"`.

- gene_sep:

  Delimiter separating genes in `gene_column` (default `"/"`).

- min_pathways:

  Minimum number of pathways a gene must appear in to be shown. Set to
  `2` to display only hub genes shared across terms. Default `1` (show
  all genes).

- max_genes:

  Maximum number of genes to display (default `50`). A safety cap for
  readability.

- gene_order_by:

  Order of gene sectors in the chord plot: `"none"` (default),
  `"foldchange"` (descending log2FC), or `"abs_foldchange"` (descending
  absolute log2FC). Fold-change based ordering requires `vista_obj` +
  `comparison`.

- gene_id_column:

  Column in `rowData(vista_obj)` used to map enrichment gene IDs to
  `vista_obj` rownames (for FC lookup).

- display_id:

  Column in `rowData(vista_obj)` providing display-friendly gene names.

- color_by:

  How to colour chords: `"foldchange"` (continuous gradient),
  `"regulation"` (Up / Down / Other), or `"pathway"` (source pathway).
  Falls back to `"pathway"` when `vista_obj` is `NULL`.

- up_color:

  Colour for up-regulated genes (default `"#D73027"`).

- down_color:

  Colour for down-regulated genes (default `"#1A9850"`).

- other_color:

  Colour for non-significant genes (default `"grey70"`).

- pathway_colors:

  Optional named vector of colours for pathway sectors. When `NULL`,
  colours are generated from an HCL palette.

- transparency:

  Chord transparency, 0–1 (default `0.4`).

- gap_degree:

  Gap between sectors in degrees (default `2`).

- label_cex:

  Text size for sector labels (default `0.7`).

- title:

  Optional plot title.

## Value

Invisibly returns a list with:

- gene_data:

  Tibble of genes with pathway membership and (optionally) fold-change
  values.

- hub_genes:

  Character vector of genes appearing in two or more pathways.

The chord diagram is drawn as a side effect.

## Details

The plot reveals **hub genes** (those driving multiple enriched terms)
and pathway redundancy (terms sharing many genes). This complements
[`get_enrichment_plot()`](get_enrichment_plot.md) (which shows
significance) and [`get_pathway_heatmap()`](get_pathway_heatmap.md)
(which shows expression patterns).

## Examples

``` r
if (FALSE) { # \dontrun{
msig <- get_msigdb_enrichment(
  vista, sample_comparison = "treatment_VS_control",
  regulation = "Up", from_type = "ENSEMBL"
)

# Simple: pathway-coloured chords
get_enrichment_chord(msig)

# With fold-change colouring
get_enrichment_chord(
  msig, vista_obj = vista,
  comparison = "treatment_VS_control",
  color_by = "foldchange"
)

# Hub genes only (shared across 2+ pathways)
get_enrichment_chord(msig, min_pathways = 2)
} # }
```
