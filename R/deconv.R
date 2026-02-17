.infer_gene_id_type <- function(gene_ids) {
  ids <- as.character(gene_ids)
  ids <- ids[!is.na(ids) & nzchar(ids)]
  if (!length(ids)) {
    return("symbol")
  }

  frac_ensembl_symbol <- mean(grepl("^ENS[A-Z]*G[0-9]+:.*", ids))
  frac_ensembl <- mean(grepl("^ENS[A-Z]*G[0-9]+(\\.[0-9]+)?$", ids))

  if (frac_ensembl_symbol >= 0.5) {
    return("ensembl_symbol")
  }
  if (frac_ensembl >= 0.5) {
    return("ensembl")
  }
  "symbol"
}

.load_xcell2_reference <- function(reference = NULL) {
  if (!is.null(reference) && !is.character(reference)) {
    return(reference)
  }

  candidates <- c(
    reference,
    "BlueprintEncode.xCell2Ref",
    "LM22.xCell2Ref",
    "TabulaSapiensBlood.xCell2Ref",
    "TabulaSapiensAll.xCell2Ref",
    "PanCancerImmune.xCell2Ref",
    "TMECompendium.xCell2Ref",
    "Immgen.xCell2Ref",
    "MouseRNAseqData.xCell2Ref",
    "DICE_demo.xCell2Ref"
  )
  candidates <- unique(candidates[!is.na(candidates) & nzchar(candidates)])

  ref_env <- new.env(parent = emptyenv())
  for (obj_name in candidates) {
    try(utils::data(list = obj_name, package = "xCell2", envir = ref_env), silent = TRUE)
    if (exists(obj_name, envir = ref_env, inherits = FALSE)) {
      return(get(obj_name, envir = ref_env, inherits = FALSE))
    }
  }

  NULL
}

.normalize_xcell2_scores <- function(scores, sample_names) {
  df <- as.data.frame(scores)
  sample_names <- as.character(sample_names)

  # Keep rows as samples when possible.
  if (!is.null(colnames(df)) && all(sample_names %in% colnames(df)) && !all(sample_names %in% rownames(df))) {
    df <- as.data.frame(t(as.matrix(df)))
  }

  if (!is.null(rownames(df)) && all(sample_names %in% rownames(df))) {
    df <- df[sample_names, , drop = FALSE]
  } else if (nrow(df) == length(sample_names)) {
    rownames(df) <- sample_names
  } else if (!is.null(colnames(df)) && all(sample_names %in% colnames(df))) {
    df <- as.data.frame(t(as.matrix(df)))
    rownames(df) <- sample_names
  } else {
    cli::cli_warn(
      c(
        "Could not confidently align deconvolution scores to samples.",
        "i" = "Expected {.val {length(sample_names)}} samples; got {.val {nrow(df)}} rows."
      )
    )
  }

  df
}

.collapse_ensembl_symbol_ids <- function(mat) {
  mat %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::mutate(gene_symbol = sub("^.+:(.+)$", "\\1", gene)) %>%
    dplyr::group_by(gene_symbol) %>%
    dplyr::summarise(dplyr::across(dplyr::where(is.numeric), mean), .groups = "drop") %>%
    tibble::column_to_rownames("gene_symbol") %>%
    as.matrix()
}

.collapse_rowdata_symbols <- function(x, mat) {
  rd <- tryCatch(SummarizedExperiment::rowData(x), error = function(e) NULL)
  if (is.null(rd) || !"SYMBOL" %in% colnames(rd)) {
    return(NULL)
  }

  symbols <- as.character(rd$SYMBOL)
  names(symbols) <- rownames(rd)
  keep <- rownames(mat) %in% names(symbols)
  if (!any(keep)) {
    return(NULL)
  }

  out <- as.data.frame(mat[keep, , drop = FALSE]) %>%
    tibble::rownames_to_column("gene_id") %>%
    dplyr::mutate(symbol = symbols[gene_id]) %>%
    dplyr::filter(!is.na(symbol), nzchar(symbol)) %>%
    dplyr::group_by(symbol) %>%
    dplyr::summarise(dplyr::across(dplyr::where(is.numeric), mean), .groups = "drop")

  if (!nrow(out)) {
    return(NULL)
  }

  out %>%
    tibble::column_to_rownames("symbol") %>%
    as.matrix()
}

#' Run Cell Deconvolution on Bulk RNA-seq from VISTA Object
#'
#' Estimates cell-type proportions in bulk RNA-seq using single-cell reference or xCell2.
#'
#' @param x A VISTA object.
#' @param method Deconvolution method. Currently only `"xCell2"` is supported.
#' @param single_cell_reference Reserved for future reference-based methods (ignored).
#' @param reference_labels Reserved for future reference-based methods (ignored).
#' @param gene_id_type Type of gene identifiers: `"auto"`, `"symbol"`, `"ensembl"`, or `"ensembl_symbol"`.
#' @param xcell2_reference Optional xCell2 reference object or dataset name
#'   (e.g., `"DICE_demo.xCell2Ref"`). Used when xCell2 exposes `xCell2Analysis()`.
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
    gene_id_type = c("auto", "symbol", "ensembl", "ensembl_symbol"),
    xcell2_reference = NULL,
    transform = c("log2", "raw"),
    ...
) {
  stopifnot(inherits(x, "VISTA"))

  method <- match.arg(method)
  gene_id_type_supplied <- !missing(gene_id_type)
  gene_id_type <- match.arg(gene_id_type)
  transform <- match.arg(transform)

  bulk_expr <- SummarizedExperiment::assay(x, "norm_counts")
  if (!gene_id_type_supplied || identical(gene_id_type, "auto")) {
    gene_id_type <- .infer_gene_id_type(rownames(bulk_expr))
  }

  if (method == "xCell2") {
    if (!requireNamespace("xCell2", quietly = TRUE)) {
      cli::cli_abort("Package {.pkg xCell2} must be installed to run `method = 'xCell2'`.")
    }
    # If gene_id_type is ensembl_symbol, extract gene symbols only
    if (gene_id_type == "ensembl_symbol") {
      bulk_expr <- .collapse_ensembl_symbol_ids(bulk_expr)
      gene_id_type <- "symbol"
    }

    xcell2_ns <- asNamespace("xCell2")
    xcell2_exports <- getNamespaceExports("xCell2")
    has_score_fun <- "xCell2Score" %in% xcell2_exports ||
      exists("xCell2Score", envir = xcell2_ns, inherits = FALSE, mode = "function")
    has_analysis_fun <- "xCell2Analysis" %in% xcell2_exports ||
      exists("xCell2Analysis", envir = xcell2_ns, inherits = FALSE, mode = "function")

    if (has_score_fun) {
      score_fun <- getFromNamespace("xCell2Score", "xCell2")
      scores <- tryCatch(
        score_fun(expr = bulk_expr, gene_id_type = gene_id_type, ...),
        error = function(e) e
      )
      if (inherits(scores, "error")) {
        cli::cli_abort(
          c(
            "xCell2 scoring failed with `xCell2Score()`.",
            "x" = conditionMessage(scores)
          )
        )
      }
    } else if (has_analysis_fun) {
      analysis_fun <- getFromNamespace("xCell2Analysis", "xCell2")
      xcell2_obj <- .load_xcell2_reference(xcell2_reference)
      if (is.null(xcell2_obj)) {
        cli::cli_abort(
          c(
            "Could not locate a default xCell2 reference object.",
            "i" = "Provide one explicitly via {.arg xcell2_reference}."
          )
        )
      }

      dots <- list(...)
      dots$xcell2object <- dots$xcell2object %||% xcell2_obj
      scores <- tryCatch(
        do.call(analysis_fun, c(list(mixture = bulk_expr), dots)),
        error = function(e) e
      )

      # Fallback: if IDs are not symbols, retry using rowData SYMBOL mapping when available.
      if (inherits(scores, "error") && gene_id_type != "symbol") {
        bulk_symbol <- .collapse_rowdata_symbols(x, bulk_expr)
        if (!is.null(bulk_symbol)) {
          scores <- tryCatch(
            do.call(analysis_fun, c(list(mixture = bulk_symbol), dots)),
            error = function(e) e
          )
        }
      }

      if (inherits(scores, "error")) {
        cli::cli_abort(
          c(
            "xCell2 scoring failed with `xCell2Analysis()`.",
            "x" = conditionMessage(scores),
            "i" = "If needed, provide a compatible xCell2 reference via {.arg xcell2_reference}."
          )
        )
      }
    } else {
      cli::cli_abort(
        c(
          "Unsupported xCell2 API detected.",
          "i" = "Expected `xCell2Score()` or `xCell2Analysis()` in package {.pkg xCell2}."
        )
      )
    }

    meta <- S4Vectors::metadata(x)
    meta$cell_fractions <- .normalize_xcell2_scores(scores, sample_names = colnames(bulk_expr))
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
