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

counts_matrix <- SummarizedExperiment::assay(airway, "counts")
count_data <- as.data.frame(counts_matrix, check.names = FALSE, stringsAsFactors = FALSE)
count_data <- tibble::rownames_to_column(count_data, var = "gene_id")
count_data <- tibble::as_tibble(count_data)

save(count_data, file = "data/count_data.rda", compress = "bzip2")
