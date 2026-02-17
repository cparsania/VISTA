if (!requireNamespace("airway", quietly = TRUE)) {
  stop("Package 'airway' is required to generate example data.")
}
if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
  stop("Package 'SummarizedExperiment' is required to generate example data.")
}
if (!requireNamespace("tibble", quietly = TRUE)) {
  stop("Package 'tibble' is required to generate example data.")
}

data("airway", package = "airway")

sample_metadata <- as.data.frame(
  SummarizedExperiment::colData(airway),
  stringsAsFactors = FALSE
)

# Keep metadata columns as plain character where applicable.
is_factor_col <- vapply(sample_metadata, is.factor, logical(1))
sample_metadata[is_factor_col] <- lapply(sample_metadata[is_factor_col], as.character)

sample_metadata$sample_names <- as.character(sample_metadata$Run)

# Keep grouping aligned with original airway metadata: dex {untrt,trt}
sample_metadata$cond_long <- ifelse(sample_metadata$dex == "trt", "treatment1", "control")
sample_metadata$cond_short <- ifelse(sample_metadata$cond_long == "control", "CTRL", "TREAT1")
sample_metadata$groups <- sample_metadata$cond_long
sample_metadata <- tibble::as_tibble(sample_metadata)

save(sample_metadata, file = "data/sample_metadata.rda", compress = "bzip2")
