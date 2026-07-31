# Display --------------------------------------------------------------------
#
# S4 objects display through show(). VISTA defined no show() method, so it fell
# back to SummarizedExperiment's, and every VISTA object introduced itself as
# "class: SummarizedExperiment" with none of its analysis state visible.

#' @keywords internal
#' @noRd
.vista_show_summary_lines <- function(object) {
  md <- S4Vectors::metadata(object)

  group_line <- NULL
  group_info <- md$group
  if (is.list(group_info) && !is.null(group_info$column)) {
    cd <- SummarizedExperiment::colData(object)
    levs <- if (group_info$column %in% colnames(cd)) {
      unique(as.character(cd[[group_info$column]]))
    } else {
      character()
    }
    group_line <- sprintf("%s (%s)", group_info$column, paste(levs, collapse = ", "))
  }

  comps <- tryCatch(names(comparisons(object)), error = function(e) NULL)
  comp_line <- if (length(comps)) {
    shown <- utils::head(comps, 3)
    more <- length(comps) - length(shown)
    sprintf(
      "%s%s", paste(shown, collapse = ", "),
      if (more > 0) sprintf(" (+%d more)", more) else ""
    )
  } else {
    "none"
  }

  source_line <- if (!is.null(md$de_active_source)) {
    available <- md$de_available_sources
    sprintf(
      "%s%s", md$de_active_source,
      if (length(available) > 1L) sprintf(" (of %s)", paste(available, collapse = ", ")) else ""
    )
  } else {
    NULL
  }

  cuts <- md$de_cutoffs
  cutoff_line <- if (is.list(cuts) && length(cuts)) {
    sprintf(
      "|log2FC| >= %s, %s <= %s",
      cuts$log2fc %||% "?", cuts$p_value_type %||% "p", cuts$pval %||% "?"
    )
  } else {
    NULL
  }

  counts_line <- if ("counts" %in% SummarizedExperiment::assayNames(object)) {
    "available via counts()"
  } else {
    NULL
  }

  fields <- list(
    "group column" = group_line,
    "comparisons" = comp_line,
    "DE source" = source_line,
    "cutoffs" = cutoff_line,
    "raw counts" = counts_line,
    "schema" = md$vista_schema_version %||% "unknown"
  )

  keep <- vapply(
    fields,
    function(v) length(v) > 0L && nzchar(paste(v, collapse = "")),
    logical(1)
  )
  fields <- fields[keep]
  if (!length(fields)) return(character())

  sprintf("%s: %s", names(fields), vapply(fields, paste, character(1), collapse = ", "))
}

#' Display a VISTA object
#'
#' Prints the `SummarizedExperiment` header followed by a short summary of the
#' analysis state -- grouping, comparisons, active differential-expression
#' source, thresholds and schema version -- none of which was visible before.
#'
#' @param object A `VISTA` object.
#' @return `object`, invisibly. Called for its output.
#' @examples
#' v <- example_vista()
#' v
#' @name VISTA-show
#' @aliases show,VISTA-method
#' @importFrom methods show
#' @exportMethod show
setMethod("show", "VISTA", function(object) {
  se_show <- methods::selectMethod("show", "SummarizedExperiment")
  out <- utils::capture.output(se_show(methods::as(object, "SummarizedExperiment")))
  # The inherited header names the parent class; this object is a VISTA.
  out[1] <- sub("^class: SummarizedExperiment", "class: VISTA", out[1])
  cat(out, sep = "\n")

  extra <- .vista_show_summary_lines(object)
  if (length(extra)) {
    cat("-------- VISTA --------", extra, sep = "\n")
  }
  invisible(object)
})
