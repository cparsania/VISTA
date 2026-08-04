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

# TRUE when a set of labels carries no identifying information -- absent, empty,
# or the positional defaults R invents ("1", "2", ...).
#' @keywords internal
#' @noRd
.vista_labels_uninformative <- function(labels, n) {
  is.null(labels) ||
    !any(nzchar(labels)) ||
    identical(as.character(labels), as.character(seq_len(n)))
}

.normalize_xcell2_scores <- function(scores, sample_names) {
  df <- as.data.frame(scores)
  sample_names <- as.character(sample_names)

  # Keep rows as samples when possible.
  if (!is.null(colnames(df)) && all(sample_names %in% colnames(df)) && !all(sample_names %in% rownames(df))) {
    df <- as.data.frame(t(as.matrix(df)))
  }

  # Preferred: align by name, in either orientation.
  if (!is.null(rownames(df)) && all(sample_names %in% rownames(df))) {
    return(df[sample_names, , drop = FALSE])
  }
  if (!is.null(colnames(df)) && all(sample_names %in% colnames(df))) {
    df <- as.data.frame(t(as.matrix(df)))
    return(df[sample_names, , drop = FALSE])
  }

  if (nrow(df) != length(sample_names)) {
    cli::cli_abort(c(
      "Could not align deconvolution scores to samples.",
      "x" = "Expected {.val {length(sample_names)}} samples; the scores have {.val {nrow(df)}} rows.",
      "i" = "Supply an {.pkg xCell2} reference whose output is labelled with the sample identifiers."
    ))
  }

  # The row count matches but the labels do not. If the scores carry real labels
  # that simply disagree, that is contradictory evidence: assigning positionally
  # would silently attach every sample's fractions to the wrong sample. Only
  # fall back to position when there is no label information to contradict.
  if (!.vista_labels_uninformative(rownames(df), nrow(df))) {
    cli::cli_abort(c(
      "Deconvolution scores are labelled, but the labels do not match the object's samples.",
      "x" = "Scores: {.val {utils::head(rownames(df), 3)}}...",
      "i" = "Object: {.val {utils::head(sample_names, 3)}}...",
      "i" = "Refusing to align by position, which would mislabel every sample."
    ))
  }

  cli::cli_inform(
    "Deconvolution scores are unlabelled; aligning to the object's {.val {length(sample_names)}} samples by position."
  )
  rownames(df) <- sample_names
  df
}

.collapse_ensembl_symbol_ids <- function(mat) {
  # Callers pass assay(x, "norm_counts"), which is a matrix; rownames_to_column()
  # requires a data.frame.
  as.data.frame(mat) %>%
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
#' @param xcell2_min_shared_genes Optional numeric shortcut for xCell2's
#'   `minSharedGenes` argument (when supported by the installed xCell2 API).
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
    xcell2_min_shared_genes = NULL,
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

    resolve_formal <- function(fun, candidates) {
      fml <- tryCatch(names(formals(fun)), error = function(e) character())
      hit <- intersect(candidates, fml)
      if (length(hit)) hit[[1]] else NULL
    }

    if (has_score_fun) {
      score_fun <- getFromNamespace("xCell2Score", "xCell2")
      dots <- list(...)
      score_arg_names <- c("expr", "mix", "mixture", "bulk_expr", "mat", "x")
      supplied_score_mat <- any(score_arg_names %in% names(dots))
      score_mat_arg <- resolve_formal(score_fun, score_arg_names)
      if (!supplied_score_mat) {
        if (is.null(score_mat_arg)) {
          dots <- c(list(bulk_expr), dots)
        } else {
          dots[[score_mat_arg]] <- bulk_expr
        }
      }
      score_gid_arg <- resolve_formal(score_fun, c("gene_id_type", "geneIdType", "gene_type", "id_type"))
      if (!is.null(score_gid_arg) && !score_gid_arg %in% names(dots)) {
        dots[[score_gid_arg]] <- gene_id_type
      }
      if (!is.null(xcell2_min_shared_genes)) {
        min_shared_arg <- resolve_formal(score_fun, c("minSharedGenes", "min_shared_genes"))
        if (!is.null(min_shared_arg) && !min_shared_arg %in% names(dots)) {
          dots[[min_shared_arg]] <- xcell2_min_shared_genes
        }
      }

      scores <- tryCatch(
        do.call(score_fun, dots),
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
      analysis_mat_arg_names <- c("mix", "mixture", "expr", "bulk_expr", "x")
      supplied_analysis_mat <- any(analysis_mat_arg_names %in% names(dots))
      analysis_mat_arg <- resolve_formal(analysis_fun, analysis_mat_arg_names)
      if (!supplied_analysis_mat) {
        if (is.null(analysis_mat_arg)) {
          dots <- c(list(bulk_expr), dots)
        } else {
          dots[[analysis_mat_arg]] <- bulk_expr
        }
      }

      analysis_ref_arg <- resolve_formal(analysis_fun, c("xcell2object", "reference", "ref"))
      if (!is.null(analysis_ref_arg) && !analysis_ref_arg %in% names(dots)) {
        dots[[analysis_ref_arg]] <- xcell2_obj
      }

      analysis_gid_arg <- resolve_formal(analysis_fun, c("gene_id_type", "geneIdType", "gene_type", "id_type"))
      if (!is.null(analysis_gid_arg) && !analysis_gid_arg %in% names(dots)) {
        dots[[analysis_gid_arg]] <- gene_id_type
      }
      if (!is.null(xcell2_min_shared_genes)) {
        min_shared_arg <- resolve_formal(analysis_fun, c("minSharedGenes", "min_shared_genes"))
        if (!is.null(min_shared_arg) && !min_shared_arg %in% names(dots)) {
          dots[[min_shared_arg]] <- xcell2_min_shared_genes
        }
      }

      scores <- tryCatch(
        do.call(analysis_fun, dots),
        error = function(e) e
      )

      # Fallback: if IDs are not symbols, retry using rowData SYMBOL mapping when available.
      if (inherits(scores, "error") && gene_id_type != "symbol") {
        bulk_symbol <- .collapse_rowdata_symbols(x, bulk_expr)
        if (!is.null(bulk_symbol)) {
          dots_symbol <- dots
          supplied_analysis_mat <- any(analysis_mat_arg_names %in% names(dots_symbol))
          if (supplied_analysis_mat && !is.null(analysis_mat_arg) && analysis_mat_arg %in% names(dots_symbol)) {
            dots_symbol[[analysis_mat_arg]] <- bulk_symbol
          } else if (length(dots_symbol) >= 1 && is.matrix(dots_symbol[[1]]) && all(dim(dots_symbol[[1]]) == dim(bulk_expr))) {
            dots_symbol[[1]] <- bulk_symbol
          } else {
            dots_symbol <- c(list(bulk_symbol), dots_symbol)
          }
          scores <- tryCatch(
            do.call(analysis_fun, dots_symbol),
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
    validate_vista(x, level = "core", error = TRUE)
    return(x)
  }

  cli::cli_abort("Deconvolution method {.val {method}} is not yet implemented.")
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

.resolve_celltype_group_column <- function(x, group_column, colnames_df) {
  if (!is.null(group_column)) {
    if (!group_column %in% colnames_df) {
      cli::cli_abort("Column {.val {group_column}} not found in sample metadata.")
    }
    return(group_column)
  }

  gcol <- tryCatch(S4Vectors::metadata(x)$group$column, error = function(e) NULL)
  if (!is.null(gcol) && gcol %in% colnames_df) {
    return(gcol)
  }

  NULL
}

.deconv_long_table <- function(x, sample_names = NULL) {
  frac <- as.data.frame(get_cell_fractions(x), stringsAsFactors = FALSE, check.names = FALSE)

  if ("sample_names" %in% colnames(frac)) {
    sample_ids <- as.character(frac$sample_names)
    frac <- dplyr::select(frac, -sample_names)
    if (anyDuplicated(sample_ids)) {
      cli::cli_abort("Cell fractions contain duplicated {.val sample_names}; sample identifiers must be unique.")
    }
    rownames(frac) <- sample_ids
  }

  if (is.null(rownames(frac)) || any(!nzchar(rownames(frac)))) {
    sample_ids <- colnames(SummarizedExperiment::assay(x, "norm_counts"))
    if (nrow(frac) == length(sample_ids)) {
      # No labels to align by, so position is the only option available.
      cli::cli_inform(
        "Cell fractions carry no sample labels; aligning to the object's samples by position."
      )
      rownames(frac) <- sample_ids
    } else {
      cli::cli_abort("Cell fractions must include sample rownames (or a {.val sample_names} column).")
    }
  }

  numeric_cols <- vapply(frac, is.numeric, logical(1))
  if (!any(numeric_cols)) {
    cli::cli_abort("Cell fractions do not contain numeric cell-type score columns.")
  }
  dropped_cols <- names(frac)[!numeric_cols]
  if (length(dropped_cols) > 0) {
    cli::cli_warn(c(
      "Dropping non-numeric columns from cell fractions.",
      "i" = "Dropped: {.val {dropped_cols}}"
    ))
    frac <- frac[, numeric_cols, drop = FALSE]
  }

  if (!is.null(sample_names)) {
    sample_names <- as.character(sample_names)
    keep <- intersect(sample_names, rownames(frac))
    if (!length(keep)) {
      cli::cli_abort("None of the requested {.arg sample_names} are present in cell fractions.")
    }
    frac <- frac[keep, , drop = FALSE]
  }

  long <- frac %>%
    tibble::rownames_to_column("sample") %>%
    tidyr::pivot_longer(-sample, names_to = "cell_type", values_to = "score")

  meta <- as.data.frame(SummarizedExperiment::colData(x))
  meta$.sample_id <- rownames(meta)

  dplyr::left_join(long, meta, by = c("sample" = ".sample_id"))
}

.choose_celltypes <- function(df, cell_types = NULL, top_n = NULL) {
  available <- unique(as.character(df$cell_type))
  if (!is.null(cell_types)) {
    missing <- setdiff(cell_types, available)
    if (length(missing) > 0) {
      cli::cli_abort("Requested cell types not found: {.val {missing}}")
    }
    return(unique(as.character(cell_types)))
  }

  if (!is.null(top_n)) {
    if (!is.numeric(top_n) || length(top_n) != 1L || top_n <= 0) {
      cli::cli_abort("{.arg top_n} must be a positive number when provided.")
    }
    top_n <- as.integer(top_n)
    return(
      df %>%
        dplyr::group_by(cell_type) %>%
        dplyr::summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop") %>%
        dplyr::arrange(dplyr::desc(mean_score)) %>%
        dplyr::slice_head(n = top_n) %>%
        dplyr::pull(cell_type)
    )
  }

  available
}

#' Plot cell-type composition as stacked bars
#'
#' @param x A VISTA object.
#' @param group_column Optional column in `sample_info(x)` used to facet/order samples.
#'   If `NULL`, uses the active VISTA group column when available.
#' @param sample_names Optional character vector of sample names to include.
#' @param base_size Base font size.
#' @param cell_types Optional character vector of cell types to keep.
#' @param top_n Optional top-N cell types by mean score (ignored when `cell_types` is provided).
#' @param collapse_other Logical; collapse non-selected cell types into `"Other"`.
#' @param normalize One of `"sample"` (default; per-sample relative scores) or `"none"`.
#' @param facet_by Faceting mode: `"group"` (default) or `"none"`.
#'
#' @return A ggplot object.
#' @examples
#' mat <- matrix(rpois(20, lambda = 20), nrow = 5)
#' rownames(mat) <- paste0("gene", seq_len(5))
#' colnames(mat) <- paste0("sample", seq_len(4))
#' se <- SummarizedExperiment::SummarizedExperiment(
#'   assays = list(norm_counts = mat),
#'   colData = S4Vectors::DataFrame(
#'     cond = c("A", "A", "B", "B"),
#'     row.names = colnames(mat)
#'   ),
#'   rowData = S4Vectors::DataFrame(
#'     gene_id = rownames(mat),
#'     row.names = rownames(mat)
#'   )
#' )
#' v <- as_vista(se, group_column = "cond")
#' md <- S4Vectors::metadata(v)
#' md$cell_fractions <- data.frame(
#'   fibroblast = c(0.2, 0.3, 0.4, 0.5),
#'   epithelial = c(0.8, 0.7, 0.6, 0.5),
#'   row.names = colnames(mat)
#' )
#' S4Vectors::metadata(v) <- md
#' get_celltype_barplot(v, group_column = "cond")
#' @export
get_celltype_barplot <- function(x,
                                 group_column = NULL,
                                 sample_names = NULL,
                                 base_size = 12,
                                 cell_types = NULL,
                                 top_n = NULL,
                                 collapse_other = TRUE,
                                 normalize = c("sample", "none"),
                                 facet_by = c("group", "none")) {
  stopifnot(inherits(x, "VISTA"))
  normalize <- match.arg(normalize)
  facet_by <- match.arg(facet_by)

  df <- .deconv_long_table(x, sample_names = sample_names)
  gcol <- .resolve_celltype_group_column(x, group_column, colnames(df))
  df$group <- if (!is.null(gcol)) as.character(df[[gcol]]) else "All"

  selected <- .choose_celltypes(df, cell_types = cell_types, top_n = top_n)

  if (isTRUE(collapse_other) && length(selected) < length(unique(df$cell_type))) {
    df <- df %>%
      dplyr::mutate(cell_type = ifelse(cell_type %in% selected, as.character(cell_type), "Other")) %>%
      dplyr::group_by(sample, group, cell_type) %>%
      dplyr::summarise(score = sum(score, na.rm = TRUE), .groups = "drop")
  } else {
    df <- df %>% dplyr::filter(cell_type %in% selected)
  }

  if (normalize == "sample") {
    df <- df %>%
      dplyr::group_by(sample) %>%
      dplyr::mutate(
        score = {
          denom <- sum(score, na.rm = TRUE)
          if (is.finite(denom) && denom > 0) score / denom else 0
        }
      ) %>%
      dplyr::ungroup()
  }

  sample_levels <- df %>%
    dplyr::distinct(sample, group) %>%
    dplyr::arrange(group, sample) %>%
    dplyr::pull(sample)
  df$sample <- factor(df$sample, levels = sample_levels)

  cell_levels <- df %>%
    dplyr::group_by(cell_type) %>%
    dplyr::summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(mean_score)) %>%
    dplyr::pull(cell_type)
  if ("Other" %in% cell_levels) {
    cell_levels <- c(setdiff(cell_levels, "Other"), "Other")
  }
  df$cell_type <- factor(df$cell_type, levels = cell_levels)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = sample, y = score, fill = cell_type)) +
    ggplot2::geom_col(width = 0.9) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
    ggplot2::xlab("Sample") +
    ggplot2::ylab(ifelse(normalize == "sample", "Relative deconvolution score", "Deconvolution score")) +
    ggplot2::labs(fill = "Cell type")

  if (facet_by == "group" && length(unique(df$group)) > 1) {
    p <- p + ggplot2::facet_wrap(~group, scales = "free_x")
  }

  p
}

#' Plot group-level deconvolution scores as dot plot
#'
#' @param x A VISTA object.
#' @param group_column Column in `sample_info(x)` that defines groups.
#'   If `NULL`, uses the active VISTA group column when available.
#' @param cell_types Optional character vector of cell types to include.
#' @param top_n Number of top cell types by mean score when `cell_types` is `NULL`.
#' @param summary_fun One of `"mean"` or `"median"` for group summary.
#' @param errorbar Error-bar type: `"se"`, `"sd"`, or `"none"`.
#' @param error Deprecated; use `errorbar`.
#' @param add_points Logical; overlay sample-level jittered points.
#' @param point_size Point size for summary points.
#' @param base_size Base font size.
#'
#' @return A ggplot object.
#' @examples
#' mat <- matrix(rpois(20, lambda = 20), nrow = 5)
#' rownames(mat) <- paste0("gene", seq_len(5))
#' colnames(mat) <- paste0("sample", seq_len(4))
#' se <- SummarizedExperiment::SummarizedExperiment(
#'   assays = list(norm_counts = mat),
#'   colData = S4Vectors::DataFrame(
#'     cond = c("A", "A", "B", "B"),
#'     row.names = colnames(mat)
#'   ),
#'   rowData = S4Vectors::DataFrame(
#'     gene_id = rownames(mat),
#'     row.names = rownames(mat)
#'   )
#' )
#' v <- as_vista(se, group_column = "cond")
#' md <- S4Vectors::metadata(v)
#' md$cell_fractions <- data.frame(
#'   fibroblast = c(0.2, 0.3, 0.4, 0.5),
#'   epithelial = c(0.8, 0.7, 0.6, 0.5),
#'   row.names = colnames(mat)
#' )
#' S4Vectors::metadata(v) <- md
#' get_celltype_group_dotplot(v, group_column = "cond")
#' @export
get_celltype_group_dotplot <- function(x,
                                       group_column = NULL,
                                       cell_types = NULL,
                                       top_n = 12,
                                       summary_fun = c("mean", "median"),
                                       errorbar = c("se", "sd", "none"),
                                       add_points = TRUE,
                                       point_size = 2.5,
                                       base_size = 12,
                                       error = NULL) {
  stopifnot(inherits(x, "VISTA"))
  summary_fun <- match.arg(summary_fun)
  errorbar <- match.arg(errorbar)
  if (!is.null(error)) {
    errorbar <- .vista_deprecate_arg(
      old = "error", new = "errorbar",
      value = match.arg(error, c("se", "sd", "none")),
      fun = "get_celltype_group_dotplot"
    )
  }
  error <- errorbar

  df <- .deconv_long_table(x)
  gcol <- .resolve_celltype_group_column(x, group_column, colnames(df))
  if (is.null(gcol)) {
    cli::cli_abort("Please provide {.arg group_column}; no default group column is available.")
  }
  df$group <- as.character(df[[gcol]])

  selected <- .choose_celltypes(df, cell_types = cell_types, top_n = top_n)
  df <- df %>% dplyr::filter(cell_type %in% selected)

  sum_tbl <- df %>%
    dplyr::group_by(group, cell_type) %>%
    dplyr::summarise(
      n = dplyr::n(),
      mean_score = mean(score, na.rm = TRUE),
      median_score = stats::median(score, na.rm = TRUE),
      sd_score = stats::sd(score, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(sd_score = ifelse(is.na(sd_score), 0, sd_score))

  sum_tbl$value <- if (summary_fun == "mean") {
    sum_tbl$mean_score
  } else {
    sum_tbl$median_score
  }

  sum_tbl$err <- if (error == "se") {
    sum_tbl$sd_score / sqrt(pmax(sum_tbl$n, 1))
  } else if (error == "sd") {
    sum_tbl$sd_score
  } else {
    0
  }

  cell_levels <- sum_tbl %>%
    dplyr::group_by(cell_type) %>%
    dplyr::summarise(m = mean(value, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(m)) %>%
    dplyr::pull(cell_type)
  sum_tbl$cell_type <- factor(sum_tbl$cell_type, levels = cell_levels)
  df$cell_type <- factor(df$cell_type, levels = cell_levels)

  dodge <- ggplot2::position_dodge(width = 0.6)
  p <- ggplot2::ggplot(sum_tbl, ggplot2::aes(x = cell_type, y = value, color = group, group = group)) +
    ggplot2::geom_point(position = dodge, size = point_size)

  if (error != "none") {
    p <- p + ggplot2::geom_errorbar(
      ggplot2::aes(ymin = value - err, ymax = value + err),
      width = 0.2,
      position = dodge
    )
  }

  if (isTRUE(add_points)) {
    p <- p + ggplot2::geom_point(
      data = df,
      ggplot2::aes(x = cell_type, y = score, color = group),
      position = ggplot2::position_jitterdodge(jitter.width = 0.15, dodge.width = 0.6),
      alpha = 0.35,
      size = point_size * 0.65,
      inherit.aes = FALSE
    )
  }

  cols <- tryCatch(group_colors(x), error = function(e) NULL)
  if (!is.null(cols) && length(cols) > 0) {
    keep <- intersect(names(cols), unique(sum_tbl$group))
    if (length(keep) > 0) {
      p <- p + ggplot2::scale_color_manual(values = cols[keep], drop = FALSE)
    }
  }

  p +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::labs(
      x = "Cell type",
      y = paste0(summary_fun, " deconvolution score"),
      color = "Group"
    )
}

#' Plot cell-type deconvolution heatmap
#'
#' @param x A VISTA object.
#' @param group_column Optional grouping column from `sample_info(x)`.
#'   If provided and `cluster_columns = FALSE`, samples are ordered by this group.
#' @param sample_names Optional character vector of sample names to include.
#' @param cell_types Optional character vector of cell types to include.
#' @param top_n Number of top cell types by mean score when `cell_types` is `NULL`.
#' @param transform One of `"none"`, `"zscore"`, or `"log1p"`.
#' @param cluster_rows Logical; hierarchical cluster cell types.
#' @param cluster_columns Logical; hierarchical cluster samples.
#' @param label Logical; overlay numeric values on tiles.
#' @param base_size Base font size.
#' @param return_type One of `"plot"` (default), `"data"`, or `"both"`. The
#'   legacy value `"matrix"` is still accepted and warns.
#'
#' @return A ggplot object, matrix, or list depending on `return_type`.
#' @examples
#' mat <- matrix(rpois(20, lambda = 20), nrow = 5)
#' rownames(mat) <- paste0("gene", seq_len(5))
#' colnames(mat) <- paste0("sample", seq_len(4))
#' se <- SummarizedExperiment::SummarizedExperiment(
#'   assays = list(norm_counts = mat),
#'   colData = S4Vectors::DataFrame(
#'     cond = c("A", "A", "B", "B"),
#'     row.names = colnames(mat)
#'   ),
#'   rowData = S4Vectors::DataFrame(
#'     gene_id = rownames(mat),
#'     row.names = rownames(mat)
#'   )
#' )
#' v <- as_vista(se, group_column = "cond")
#' md <- S4Vectors::metadata(v)
#' md$cell_fractions <- data.frame(
#'   fibroblast = c(0.2, 0.3, 0.4, 0.5),
#'   epithelial = c(0.8, 0.7, 0.6, 0.5),
#'   row.names = colnames(mat)
#' )
#' S4Vectors::metadata(v) <- md
#' get_celltype_heatmap(v, group_column = "cond")
#' @export
get_celltype_heatmap <- function(x,
                                 group_column = NULL,
                                 sample_names = NULL,
                                 cell_types = NULL,
                                 top_n = 20,
                                 transform = c("none", "zscore", "log1p"),
                                 cluster_rows = TRUE,
                                 cluster_columns = TRUE,
                                 label = FALSE,
                                 base_size = 11,
                                 return_type = c("plot", "data", "both")) {
  stopifnot(inherits(x, "VISTA"))
  transform <- match.arg(transform)
  return_type <- .vista_resolve_return_type(
    return_type, fun = "get_celltype_heatmap", legacy = c(matrix = "data")
  )

  df <- .deconv_long_table(x, sample_names = sample_names)
  selected <- .choose_celltypes(df, cell_types = cell_types, top_n = top_n)
  df <- df %>% dplyr::filter(cell_type %in% selected)

  mat_df <- df %>%
    dplyr::select(cell_type, sample, score) %>%
    tidyr::pivot_wider(names_from = sample, values_from = score, values_fill = 0)

  mat <- as.matrix(mat_df[, -1, drop = FALSE])
  rownames(mat) <- mat_df$cell_type

  if (transform == "log1p") {
    mat <- log1p(mat)
  } else if (transform == "zscore") {
    mat <- t(scale(t(mat)))
    mat[!is.finite(mat)] <- 0
  }

  if (isTRUE(cluster_rows) && nrow(mat) > 1) {
    mat <- mat[hclust(stats::dist(mat))$order, , drop = FALSE]
  }

  if (isTRUE(cluster_columns) && ncol(mat) > 1) {
    mat <- mat[, hclust(stats::dist(t(mat)))$order, drop = FALSE]
  } else if (!isTRUE(cluster_columns)) {
    gcol <- .resolve_celltype_group_column(x, group_column, colnames(df))
    if (!is.null(gcol)) {
      ord <- df %>%
        dplyr::distinct(sample, group = .data[[gcol]]) %>%
        dplyr::arrange(group, sample) %>%
        dplyr::pull(sample)
      ord <- ord[ord %in% colnames(mat)]
      if (length(ord) > 0) mat <- mat[, ord, drop = FALSE]
    }
  }

  if (return_type == "data") {
    return(mat)
  }

  plot_df <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  colnames(plot_df) <- c("cell_type", "sample", "value")
  plot_df$sample <- factor(plot_df$sample, levels = colnames(mat))
  plot_df$cell_type <- factor(plot_df$cell_type, levels = rev(rownames(mat)))

  fill_title <- if (transform == "zscore") "Z-score" else if (transform == "log1p") "log1p(score)" else "Score"

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = sample, y = cell_type, fill = value)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c(name = fill_title) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )

  if (isTRUE(label)) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.2f", value)),
      size = max(2.5, base_size / 4)
    )
  }

  if (return_type == "both") {
    return(list(plot = p, matrix = mat))
  }

  p
}
