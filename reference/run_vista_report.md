# Generate a publication-ready VISTA workflow report

Builds a `VISTA` object (or uses a precomputed one), computes a
comprehensive single-comparison analysis panel, exports publication
assets, and renders an automated Quarto HTML/PDF report from YAML-driven
parameters.

## Usage

``` r
run_vista_report(config, output_file = "vista-report.html")
```

## Arguments

- config:

  Path to a YAML config file (see template at
  `inst/reports/vista-report-template.yml`) or a named list of
  parameters.

- output_file:

  Optional output filename overriding `config$output_file`.

## Value

Invisibly, the normalized output report path.

## Details

The report focuses on one differential comparison (`primary_comparison`)
and includes:

- QC plots (PCA, MDS, correlation heatmap),

- DE plots (volcano, MA, DEG bar/pie/donut summaries),

- expression-focused views (boxplot, fold-change barplot, expression
  heatmap),

- enrichment outputs (MSigDB/GO/KEGG when available),

- downloadable artifacts (tables + plots + optional zip bundle),

- interactive HTML tables via DT when installed.

## Examples

``` r
if (FALSE) { # \dontrun{
data('count_data', package = 'VISTA')
data('sample_metadata', package = 'VISTA')
cfg <- list(
  counts = count_data[1:100, ],
  sample_info = sample_metadata[1:6, ],
  column_geneid = 'gene_id',
  group_column = 'cond_long',
  group_numerator = 'treatment1',
  group_denominator = 'control',
  include_msigdb = FALSE, include_go = FALSE, include_kegg = FALSE
)
if (requireNamespace('quarto', quietly = TRUE)) {
  out <- tempfile(fileext = '.html')
  try(run_vista_report(cfg, output_file = out), silent = TRUE)
}
} # }
```
