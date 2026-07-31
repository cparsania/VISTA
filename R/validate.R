#' Validate a VISTA object
#'
#' Performs structural and metadata contract checks on a `VISTA` object.
#' Use `level = "core"` for invariants that all VISTA objects must satisfy,
#' or `level = "full"` to additionally validate method-specific metadata used
#' by advanced workflows.
#'
#' @param x A `VISTA` object.
#' @param level Validation depth. One of `"core"` or `"full"`.
#' @param error Logical; if `TRUE`, aborts on validation failures. If `FALSE`,
#'   returns the validation report and emits a warning when issues are found.
#'
#' @return Invisibly returns a list with fields `valid`, `level`, and `issues`.
#' @examples
#' mat <- matrix(rnorm(60), nrow = 10)
#' rownames(mat) <- paste0("gene", seq_len(nrow(mat)))
#' colnames(mat) <- paste0("sample", seq_len(ncol(mat)))
#' se <- SummarizedExperiment::SummarizedExperiment(
#'   assays = list(norm_counts = mat),
#'   colData = S4Vectors::DataFrame(
#'     cond = rep(c("A", "B"), each = 3),
#'     row.names = colnames(mat)
#'   ),
#'   rowData = S4Vectors::DataFrame(
#'     gene_id = rownames(mat),
#'     row.names = rownames(mat)
#'   )
#' )
#' v <- as_vista(se, group_column = "cond")
#' validate_vista(v)
#' @export
validate_vista <- function(x, level = c("core", "full"), error = TRUE) {
  level <- match.arg(level)

  issues <- .validate_vista_core(x)
  if (identical(level, "full") && inherits(x, "VISTA")) {
    issues <- c(issues, .validate_vista_full(x))
  }

  out <- list(
    valid = length(issues) == 0,
    level = level,
    issues = unique(issues)
  )

  if (!out$valid) {
    bullets <- stats::setNames(out$issues, rep("x", length(out$issues)))
    if (isTRUE(error)) {
      cli::cli_abort(c("VISTA validation failed.", bullets))
    } else {
      cli::cli_warn(c("VISTA validation reported issues.", bullets))
    }
  }

  invisible(out)
}

#' @keywords internal
.validate_vista_core <- function(x) {
  issues <- character()

  if (!inherits(x, "VISTA")) {
    return("Object must inherit from class 'VISTA'.")
  }

  valid_msg <- methods::validObject(x, test = TRUE)
  if (!isTRUE(valid_msg)) {
    issues <- c(issues, as.character(valid_msg))
  }

  mat <- tryCatch(SummarizedExperiment::assay(x, "norm_counts"), error = function(e) NULL)
  if (is.null(mat)) {
    issues <- c(issues, "Assay 'norm_counts' is missing.")
    return(unique(issues))
  }

  si <- as.data.frame(SummarizedExperiment::colData(x), stringsAsFactors = FALSE)
  rd <- as.data.frame(SummarizedExperiment::rowData(x), stringsAsFactors = FALSE)
  md <- S4Vectors::metadata(x)

  if (ncol(mat) != nrow(si)) {
    issues <- c(issues, "ncol(norm_counts) must equal nrow(colData).")
  }
  if (nrow(mat) != nrow(rd)) {
    issues <- c(issues, "nrow(norm_counts) must equal nrow(rowData).")
  }

  if ("sample_names" %in% colnames(si) && !identical(as.character(si$sample_names), colnames(mat))) {
    issues <- c(issues, "sample_info$sample_names must match colnames(norm_counts).")
  }

  if ("gene_id" %in% colnames(rd) && !identical(as.character(rd$gene_id), rownames(mat))) {
    issues <- c(issues, "row_data$gene_id should match rownames(norm_counts).")
  }

  if (!is.character(md$vista_schema_version) || length(md$vista_schema_version) != 1L || !nzchar(md$vista_schema_version)) {
    issues <- c(issues, "metadata(x)$vista_schema_version must be a non-empty character scalar.")
  } else {
    schema_state <- .vista_schema_compare(x)
    current_schema <- .VISTA_SCHEMA_VERSION
    stored_schema <- md$vista_schema_version
    if (identical(schema_state, "newer")) {
      # A newer layout may carry semantics this version would misread, so this
      # is an issue rather than a note.
      issues <- c(issues, sprintf(
        "metadata(x)$vista_schema_version ('%s') is newer than this version of VISTA supports ('%s'); upgrade the package.",
        md$vista_schema_version, .VISTA_SCHEMA_VERSION
      ))
    } else if (identical(schema_state, "older")) {
      cli::cli_inform(c(
        "This VISTA object uses schema {.val {stored_schema}}; the current schema is {.val {current_schema}}.",
        "i" = "Run {.code updateObject(x)} to migrate it."
      ))
    }
  }

  if (is.list(md$de_cutoffs)) {
    if (!is.null(md$de_cutoffs$p_value_type) && !md$de_cutoffs$p_value_type %in% c("padj", "pvalue")) {
      issues <- c(issues, "metadata(x)$de_cutoffs$p_value_type must be one of: padj, pvalue.")
    }
  }

  # The DE tables must describe the genes actually present. This lived only in
  # .validate_vista_full(), so a stale-metadata object -- the state row
  # subsetting used to produce -- passed a core check. It stays out of
  # setValidity(), which runs on every new() and slot assignment and must remain
  # cheap and structural.
  issues <- c(issues, .validate_de_rownames(md$de_results, rownames(mat), "metadata(x)$de_results"))
  if (!is.null(md$de_results_by_method) && is.list(md$de_results_by_method)) {
    for (src in names(md$de_results_by_method)) {
      issues <- c(issues, .validate_de_rownames(
        md$de_results_by_method[[src]], rownames(mat),
        sprintf("metadata(x)$de_results_by_method[['%s']]", src)
      ))
    }
  }

  unique(issues)
}

#' @keywords internal
.validate_de_rownames <- function(tbls, rn, label) {
  if (is.null(tbls)) return(character())
  if (inherits(tbls, "SimpleList")) tbls <- as.list(tbls)
  if (!is.list(tbls)) return(character())

  out <- character()
  for (nm in names(tbls)) {
    trn <- rownames(as.data.frame(tbls[[nm]], stringsAsFactors = FALSE))
    if (is.null(trn) || !identical(trn, rn)) {
      out <- c(out, sprintf("%s[['%s']] rownames must match rownames(norm_counts).", label, nm))
    }
  }
  out
}

#' @keywords internal
.validate_vista_full <- function(x) {
  issues <- character()
  rn <- rownames(SummarizedExperiment::assay(x, "norm_counts"))
  md <- S4Vectors::metadata(x)

  as_named_list <- function(obj) {
    if (is.null(obj)) return(list())
    if (inherits(obj, "SimpleList")) obj <- as.list(obj)
    if (!is.list(obj)) return(NULL)
    obj
  }

  check_result_tables <- function(tbls, label) {
    out <- character()
    if (is.null(tbls)) return(out)
    if (inherits(tbls, "SimpleList")) tbls <- as.list(tbls)
    if (!is.list(tbls)) {
      return(sprintf("%s must be a list/SimpleList of DE tables.", label))
    }
    for (nm in names(tbls)) {
      tbl <- as.data.frame(tbls[[nm]], stringsAsFactors = FALSE)
      trn <- rownames(tbl)
      if (is.null(trn) || !identical(trn, rn)) {
        out <- c(out, sprintf("%s[['%s']] rownames must match rownames(norm_counts).", label, nm))
      }
      required_cols <- c("log2fc", "regulation")
      if (length(setdiff(required_cols, colnames(tbl))) > 0) {
        out <- c(out, sprintf("%s[['%s']] should include columns: %s.", label, nm, paste(required_cols, collapse = ", ")))
      }
    }
    out
  }

  issues <- c(issues, check_result_tables(md$de_results, "metadata(x)$de_results"))

  by_method <- as_named_list(md$de_results_by_method)
  if (!is.null(by_method) && length(by_method) > 0) {
    for (src in names(by_method)) {
      issues <- c(issues, check_result_tables(by_method[[src]], sprintf("metadata(x)$de_results_by_method[['%s']]", src)))
    }

    if (is.null(md$de_active_source) || !md$de_active_source %in% names(by_method)) {
      issues <- c(issues, "metadata(x)$de_active_source must be one of names(metadata(x)$de_results_by_method).")
    }
  }

  if (!is.null(md$comparison) && is.list(md$comparison)) {
    comp_cols <- md$comparison$colors
    active <- comparisons(x)
    active_names <- names(active)
    if (length(active_names) > 0 && !is.null(comp_cols)) {
      missing <- setdiff(active_names, names(comp_cols))
      if (length(missing) > 0) {
        issues <- c(issues, sprintf("comparison color map missing entries for: %s.", paste(missing, collapse = ", ")))
      }
    }
  }

  unique(issues)
}
