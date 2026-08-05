# Subsetting -----------------------------------------------------------------
#
# VISTA inherits SummarizedExperiment's `[`, which knows nothing about the DE
# tables in metadata(). Before 1.2.0, v[1:10, ] returned a 10-row object whose
# metadata still described all the original genes, and validObject() reported it
# as valid. Column subsetting was worse: it passed even validate_vista(level =
# "full") while leaving colour maps naming groups that no longer existed.

#' @keywords internal
#' @noRd
.vista_subset_de_tables <- function(tbls, keep) {
  if (is.null(tbls)) return(NULL)
  was_simple <- inherits(tbls, "SimpleList")
  lst <- if (was_simple) as.list(tbls) else tbls
  if (!is.list(lst)) return(tbls)

  out <- lapply(lst, function(tbl) {
    df <- as.data.frame(tbl, stringsAsFactors = FALSE)
    rn <- rownames(df)
    if (is.null(rn)) return(df)
    idx <- match(keep, rn)
    idx <- idx[!is.na(idx)]
    df[idx, , drop = FALSE]
  })

  if (was_simple) S4Vectors::SimpleList(out) else out
}

#' @keywords internal
#' @noRd
.vista_subset_summaries <- function(summaries, comps) {
  if (is.null(summaries)) return(NULL)
  was_simple <- inherits(summaries, "SimpleList")
  lst <- if (was_simple) as.list(summaries) else summaries
  if (!is.list(lst)) return(summaries)

  # Recount regulation classes from the retained DE rows so the summary keeps
  # describing the object it belongs to.
  out <- lapply(names(lst), function(nm) {
    tbl <- comps[[nm]]
    if (is.null(tbl) || !"regulation" %in% colnames(tbl)) return(lst[[nm]])
    dplyr::count(as.data.frame(tbl, stringsAsFactors = FALSE), regulation)
  })
  names(out) <- names(lst)

  if (was_simple) S4Vectors::SimpleList(out) else out
}

#' Subset a VISTA object
#'
#' @description
#' Subsets the underlying `SummarizedExperiment` and keeps the analysis metadata
#' consistent with the result.
#'
#' @details
#' Row (gene) subsetting reindexes every stored differential-expression table --
#' the active `de_results`/`de_summary` and every entry of
#' `de_results_by_method`/`de_summary_by_method` -- to the retained genes, in the
#' new row order. DEG summaries are recounted from the retained rows.
#'
#' Column (sample) subsetting leaves the DE tables untouched, because they are
#' analysis *results* rather than per-sample data; recomputing them would
#' require re-running the model. Group and comparison colour maps are pruned to
#' the groups that still have samples, and a warning names any comparison whose
#' groups are no longer represented, since plots of that comparison no longer
#' correspond to the samples in the object.
#'
#' Each subset appends an entry to `metadata(x)$provenance$subset_history`.
#'
#' @param x A `VISTA` object.
#' @param i,j Row (gene) and column (sample) subscripts.
#' @param ... Passed to the `SummarizedExperiment` method.
#' @param drop Ignored, as for `SummarizedExperiment`.
#'
#' @return A `VISTA` object.
#' @examples
#' v <- example_vista()
#'
#' # DE tables follow the genes.
#' v10 <- v[seq_len(10), ]
#' nrow(comparisons(v10)[[1]])
#'
#' # Sample subsetting prunes the group colour map.
#' v_ctrl <- v[, sample_info(v)$cond_long == "control"]
#' group_colors(v_ctrl)
#'
#' @name VISTA-subset
#' @aliases [,VISTA,ANY,ANY-method [,VISTA,ANY,ANY,ANY-method
#' @exportMethod [
setMethod("[", c("VISTA", "ANY", "ANY"), function(x, i, j, ..., drop = FALSE) {
  if (!missing(drop) && !identical(drop, FALSE)) {
    cli::cli_warn("{.arg drop} is ignored when subsetting a {.cls VISTA} object.")
  }

  md <- S4Vectors::metadata(x)
  out <- callNextMethod(x, i, j, ..., drop = FALSE)

  subset_rows <- !missing(i)
  subset_cols <- !missing(j)

  if (subset_rows) {
    keep <- rownames(out)
    md$de_results <- .vista_subset_de_tables(md$de_results, keep)
    if (!is.null(md$de_results_by_method)) {
      md$de_results_by_method <- lapply(
        md$de_results_by_method, .vista_subset_de_tables, keep = keep
      )
    }
    md$de_summary <- .vista_subset_summaries(md$de_summary, md$de_results)
    if (!is.null(md$de_summary_by_method)) {
      md$de_summary_by_method <- lapply(
        names(md$de_summary_by_method),
        function(src) {
          .vista_subset_summaries(
            md$de_summary_by_method[[src]], md$de_results_by_method[[src]]
          )
        }
      )
      names(md$de_summary_by_method) <- names(md$de_results_by_method)
    }
  }

  if (subset_cols) {
    group_info <- md$group
    if (is.list(group_info) && !is.null(group_info$column)) {
      cd <- SummarizedExperiment::colData(out)
      gcol <- group_info$column
      if (gcol %in% colnames(cd)) {
        remaining <- unique(as.character(cd[[gcol]]))
        if (!is.null(group_info$colors)) {
          group_info$colors <- group_info$colors[
            intersect(names(group_info$colors), remaining)
          ]
        }
        md$group <- group_info

        # A comparison whose numerator or denominator has no samples left is
        # still stored, but its plots no longer describe this object.
        comp_names <- names(md$de_results)
        if (length(comp_names)) {
          orphaned <- comp_names[vapply(comp_names, function(nm) {
            parts <- strsplit(nm, "_VS_", fixed = TRUE)[[1]]
            length(parts) == 2L && !all(parts %in% remaining)
          }, logical(1))]
          if (length(orphaned)) {
            cli::cli_warn(c(
              "{cli::qty(orphaned)}Comparison{?s} {.val {orphaned}} {?is/are} no longer represented by the retained samples.",
              "i" = "The stored results are kept unchanged; re-run {.fn create_vista} to analyse this subset."
            ))
          }
        }
      }
    }
  }

  if (subset_rows || subset_cols) {
    prov <- md$provenance
    if (!is.list(prov)) prov <- list()
    prov$subset_history <- c(prov$subset_history, list(list(
      genes = nrow(out),
      samples = ncol(out),
      at = as.character(Sys.time())
    )))
    md$provenance <- prov
  }

  S4Vectors::metadata(out) <- md
  out
})
