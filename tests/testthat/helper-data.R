make_small_vista <- function() {
  example_vista()
}

make_benchmark_inputs <- function(n_genes = 150) {
  data("count_data", package = "VISTA", envir = environment())
  data("sample_metadata", package = "VISTA", envir = environment())

  target_groups <- c("control", "treatment1")
  sample_subset <- sample_metadata[sample_metadata$cond_long %in% target_groups, , drop = FALSE]
  counts_subset <- count_data[seq_len(n_genes), c("gene_id", sample_subset$sample_names)]

  list(
    counts = counts_subset,
    sample_info = sample_subset,
    column_geneid = "gene_id",
    group_column = "cond_long",
    group_numerator = "treatment1",
    group_denominator = "control",
    min_counts = 5,
    min_replicates = 1,
    log2fc_cutoff = 1,
    pval_cutoff = 0.05,
    p_value_type = "padj"
  )
}
