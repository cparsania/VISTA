#' @keywords internal
.build_example_vista <- function(n_genes = 150,
                                 n_per_group = 3,
                                 method = "deseq2") {
  data("count_data", package = "VISTA", envir = environment())
  data("sample_metadata", package = "VISTA", envir = environment())

  stopifnot(is.numeric(n_genes), length(n_genes) == 1L, n_genes > 0)
  stopifnot(is.numeric(n_per_group), length(n_per_group) == 1L, n_per_group > 0)

  si <- as.data.frame(sample_metadata, stringsAsFactors = FALSE)
  req_groups <- c("control", "treatment1")
  if (!all(req_groups %in% unique(si$cond_long))) {
    cli::cli_abort(
      "Built-in {.data sample_metadata} must contain groups {.val {req_groups}} in {.field cond_long}."
    )
  }

  ctrl <- si[si$cond_long == "control", , drop = FALSE]
  trt <- si[si$cond_long == "treatment1", , drop = FALSE]
  ctrl <- utils::head(ctrl, n_per_group)
  trt <- utils::head(trt, n_per_group)
  sample_subset <- rbind(ctrl, trt)

  counts_subset <- count_data[seq_len(min(n_genes, nrow(count_data))), , drop = FALSE]
  counts_subset <- counts_subset[, c("gene_id", sample_subset$sample_names), drop = FALSE]

  create_vista(
    counts = counts_subset,
    sample_info = sample_subset,
    column_geneid = "gene_id",
    group_column = "cond_long",
    group_numerator = "treatment1",
    group_denominator = "control",
    min_counts = 5,
    min_replicates = 1,
    method = method
  )
}

#' Build a small example VISTA object
#'
#' Creates a lightweight `VISTA` object from built-in package datasets
#' (`count_data` and `sample_metadata`) for use in examples and tutorials.
#' The default call returns a precomputed object to keep examples and package
#' checks fast. Non-default argument combinations fall back to rebuilding the
#' object from the packaged example inputs.
#'
#' @param n_genes Number of genes to include (default `150`).
#' @param n_per_group Number of samples per group (`control` and `treatment1`)
#'   to include (default `3`).
#' @param method Differential expression backend passed to [create_vista()].
#'
#' @return A `VISTA` object.
#' @export
#' @examples
#' v <- example_vista()
#' v
example_vista <- function(n_genes = 150,
                          n_per_group = 3,
                          method = "deseq2") {
  stopifnot(is.numeric(n_genes), length(n_genes) == 1L, n_genes > 0)
  stopifnot(is.numeric(n_per_group), length(n_per_group) == 1L, n_per_group > 0)

  if (identical(as.integer(n_genes), 150L) &&
      identical(as.integer(n_per_group), 3L) &&
      identical(method, "deseq2")) {
    return(vista_example_default)
  }

  .build_example_vista(
    n_genes = n_genes,
    n_per_group = n_per_group,
    method = method
  )
}
