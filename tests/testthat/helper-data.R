make_small_vista <- function() {
  data("count_data", package = "VISTA", envir = environment())
  data("sample_metadata", package = "VISTA", envir = environment())

  ctrl_subset <- sample_metadata[sample_metadata$cond_long == "control", ][1:3, , drop = FALSE]
  trt_subset <- sample_metadata[sample_metadata$cond_long == "treatment1", ][1:3, , drop = FALSE]
  sample_subset <- rbind(ctrl_subset, trt_subset)
  counts_subset <- count_data[seq_len(150), c("gene_id", sample_subset$sample_names)]

  create_vista(
    counts = counts_subset,
    sample_info = sample_subset,
    column_geneid = "gene_id",
    group_column = "cond_long",
    group_numerator = "treatment1",
    group_denominator = "control",
    min_counts = 5,
    min_replicates = 1
  )
}
