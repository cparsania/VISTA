#' Run Cell Deconvolution on Bulk RNA-seq from VISTA Object
#'
#' Estimates cell-type proportions in bulk RNA-seq using single-cell reference or xCell2.
#'
#' @param x A VISTA object.
#' @param method Deconvolution method. Currently only `"xCell2"` is supported.
#' @param single_cell_reference Reserved for future reference-based methods (ignored).
#' @param reference_labels Reserved for future reference-based methods (ignored).
#' @param gene_id_type Type of gene identifiers: "symbol", "ensembl", or "ensembl_symbol".
#' @param transform Expression transformation: "log2" or "raw".
#' @param ... Additional arguments passed to the specific method.
#'
#' @return VISTA object with cell_fractions added to metadata.
#' @export
run_cell_deconvolution <- function(
    x,
    method = c("xCell2"),
    single_cell_reference = NULL,
    reference_labels = NULL,
    gene_id_type = c("symbol", "ensembl", "ensembl_symbol"),
    transform = c("log2", "raw"),
    ...
) {
  stopifnot(inherits(x, "VISTA"))

  method <- match.arg(method)
  gene_id_type <- match.arg(gene_id_type)
  transform <- match.arg(transform)

  bulk_expr <- SummarizedExperiment::assay(x, "norm_counts")

  if (method == "xCell2") {
    if (!requireNamespace("xCell2", quietly = TRUE)) {
      cli::cli_abort("Package {.pkg xCell2} must be installed to run `method = 'xCell2'`.")
    }
    # If gene_id_type is ensembl_symbol, extract gene symbols only
    if (gene_id_type == "ensembl_symbol") {
      bulk_expr <- bulk_expr %>%
        tibble::rownames_to_column("gene") %>%
        dplyr::mutate(gene_symbol = sub("^.+:(.+)$", "\1", gene)) %>%
        dplyr::group_by(gene_symbol) %>%
        dplyr::summarise(dplyr::across(dplyr::where(is.numeric), mean), .groups = "drop") %>%
        tibble::column_to_rownames("gene_symbol") %>%
        as.matrix()
      gene_id_type <- "symbol"
    }
    scores <- xCell2::xCell2Score(expr = bulk_expr, gene_id_type = gene_id_type, ...)
    meta <- S4Vectors::metadata(x)
    meta$cell_fractions <- as.data.frame(scores)
    S4Vectors::metadata(x) <- meta
    return(x)
  }

  cli::cli_abort("Deconvolution method {.val {method}} is not yet implemented.")
}

#' Plot cell type proportions
#'
#' @param x A VISTA object.
#' @param group_column A column in colData(x) to group by (e.g., "condition").
#' @param samples Optional character vector of sample names to include.
#' @param font_size Base font size.
#' @return A ggplot object.
#' @export
plot_celltype_barplot <- function(x,
                                  group_column = NULL,
                                  samples = NULL,
                                  font_size = 12) {
  df <- get_cell_fractions(x) %>%
    tibble::rownames_to_column("sample") %>%
    tidyr::pivot_longer(-sample, names_to = "cell_type", values_to = "proportion")

  df <- df %>%
    dplyr::left_join(as.data.frame(SummarizedExperiment::colData(x)) %>% tibble::rownames_to_column("sample"), by = "sample")

  if (!is.null(samples)) {
    df <- df %>% dplyr::filter(sample %in% samples)
  }

  if (!is.null(group_column) && group_column %in% colnames(df)) {
    df$group <- df[[group_column]]
  } else {
    df$group <- "All"
  }

  ggplot2::ggplot(df, ggplot2::aes(x = sample, y = proportion, fill = cell_type)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::facet_wrap(~group, scales = "free_x") +
    ggplot2::theme_minimal(base_size = font_size) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
    ggplot2::ylab("Estimated Proportion") +
    ggplot2::xlab("Sample")
}

#' Retrieve stored cell fraction estimates
#'
#' @param x A VISTA object.
#' @return A data.frame of cell fraction estimates with samples in rows.
#' @export
get_cell_fractions <- function(x) {
  stopifnot(inherits(x, "VISTA"))
  meta <- S4Vectors::metadata(x)
  res <- meta$cell_fractions
  if (is.null(res)) {
    cli::cli_abort("No cell fraction estimates found in metadata(x)$cell_fractions. Run `run_cell_deconvolution()` first.")
  }
  as.data.frame(res)
}
