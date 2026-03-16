#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(VISTA)
})

data("count_data", package = "VISTA")
data("sample_metadata", package = "VISTA")

target_groups <- c("control", "treatment1")
sample_subset <- sample_metadata[sample_metadata$cond_long %in% target_groups, , drop = FALSE]
count_subset <- count_data[1:10000, c("gene_id", sample_subset$sample_names)]

report <- benchmark_vista_equivalence(
  counts = count_subset,
  sample_info = sample_subset,
  column_geneid = "gene_id",
  group_column = "cond_long",
  group_numerator = "treatment1",
  group_denominator = "control",
  methods = c("deseq2", "edger", "limma"),
  min_counts = 5,
  min_replicates = 1,
  log2fc_cutoff = 1,
  pval_cutoff = 0.05,
  p_value_type = "padj",
  tolerance = 1e-8,
  return_plots = TRUE
)

cat("\nComparison summary\n")
print(report$comparison_summary)

cat("\nVisual summary\n")
print(report$visual_summary)

if (!report$valid) {
  stop("Benchmark detected differences between VISTA and standalone backends.", call. = FALSE)
}

cat("\nAll VISTA backend checks matched the standalone pipelines.\n")
