# Deprecation infrastructure -------------------------------------------------
#
# VISTA follows the Bioconductor deprecation contract: an argument warns for one
# release, becomes defunct (errors) in the next, and is removed in the one after.
# Every timeline lives in `.vista_deprecations()` rather than at the call site,
# so a single table is the source of truth for the code, the `?VISTA-deprecated`
# man page, and the test that asserts each alias still warns.

#' @keywords internal
#' @noRd
.vista_deprecations <- function() {
  # fun, old_arg, new_arg, deprecated_in, defunct_in, removed_in, note
  #
  # INVARIANT: every row must describe an alias that already exists in the code.
  # `test-deprecation-registry.R` asserts that `old_arg` is still a formal (or
  # reaches the function through `...`) and that `new_arg` is a formal too. Add
  # a row in the same commit that adds its alias, never ahead of it.
  rows <- list(
    c("get_corr_heatmap", "show_corr_values", "label", "1.2.0", "1.4.0", "1.6.0", ""),
    c("get_corr_heatmap", "col_corr_values", "label_color", "1.2.0", "1.4.0", "1.6.0", ""),
    c("get_volcano_plot", "col_up", "colors", "1.2.0", "1.4.0", "1.6.0", "Supply colors = c(Up = ...)."),
    c("get_volcano_plot", "col_down", "colors", "1.2.0", "1.4.0", "1.6.0", "Supply colors = c(Down = ...)."),
    c("get_volcano_plot", "col_other", "colors", "1.2.0", "1.4.0", "1.6.0", "Supply colors = c(Other = ...)."),
    c("get_volcano_plot", "col_others", "colors", "1.2.0", "1.4.0", "1.6.0", "Supply colors = c(Other = ...)."),
    c("get_volcano_plot", "lab_size", "label_size", "1.2.0", "1.4.0", "1.6.0", ""),
    c("get_pca_plot", "sample.seed", "", "1.2.0", "1.4.0", "1.6.0",
      "`sample.seed` has never had any effect; PCA is deterministic."),
    c("get_pca_plot", "use_vista_colors", "use_group_colors", "1.2.0", "1.4.0", "1.6.0", ""),
    c("get_mds_plot", "use_vista_colors", "use_group_colors", "1.2.0", "1.4.0", "1.6.0", ""),
    c("get_umap_plot", "use_vista_colors", "use_group_colors", "1.2.0", "1.4.0", "1.6.0", ""),
    c("get_expression_barplot", "facet_scale", "facet_scales", "1.2.0", "1.4.0", "1.6.0", ""),
    c("get_expression_violinplot", "value_transform", "log_transform", "1.2.0", "1.4.0", "1.6.0", ""),
    c("get_expression_lineplot", "value_transform", "log_transform", "1.2.0", "1.4.0", "1.6.0", "")
  )

  out <- as.data.frame(
    do.call(rbind, rows),
    stringsAsFactors = FALSE
  )
  names(out) <- c("fun", "old_arg", "new_arg", "deprecated_in", "defunct_in", "removed_in", "note")
  out
}

#' @keywords internal
#' @noRd
.vista_deprecation_row <- function(fun, old) {
  reg <- .vista_deprecations()
  hit <- reg[reg$fun == fun & reg$old_arg == old, , drop = FALSE]
  if (nrow(hit) == 1L) hit else NULL
}

#' @keywords internal
#' @noRd
.vista_deprecate_msg <- function(fun, old, new, row) {
  target <- if (!is.null(new) && nzchar(new)) {
    sprintf("Use `%s` instead.", new)
  } else {
    "It no longer has any effect."
  }
  when <- if (!is.null(row) && nzchar(row$defunct_in)) {
    sprintf(" It becomes defunct in VISTA %s.", row$defunct_in)
  } else {
    ""
  }
  extra <- if (!is.null(row) && nzchar(row$note)) paste0(" ", row$note) else ""
  sprintf("`%s` in `%s()` is deprecated. %s%s%s", old, fun, target, when, extra)
}

#' Signal a deprecated VISTA argument
#'
#' Emits a single warning classed `c("vista_deprecated_arg",
#' "deprecatedWarning")` and returns the (optionally transformed) legacy value.
#' The `deprecatedWarning` class matches [base::.Deprecated()], so
#' `suppressWarnings()`, `tryCatch()` and Bioconductor tooling all behave as
#' expected; the `vista_deprecated_arg` subclass is what VISTA's own tests
#' assert on.
#'
#' @param old,new Character scalars naming the deprecated and replacement arguments.
#' @param value The value supplied to the deprecated argument.
#' @param fun Character scalar naming the calling function.
#' @param transform Optional function applied to `value` before returning.
#'
#' @return `transform(value)`.
#' @keywords internal
#' @noRd
.vista_deprecate_arg <- function(old, new = NULL, value, fun, transform = identity) {
  row <- .vista_deprecation_row(fun, old)
  if (is.null(new) && !is.null(row)) new <- row$new_arg

  msg <- .vista_deprecate_msg(fun, old, new, row)
  cls <- c("vista_deprecated_arg", "deprecatedWarning")

  # Throttled by default so a deprecated argument inside a loop cannot emit
  # hundreds of warnings. Tests (and users who want every occurrence) set
  # `options(vista.deprecation_frequency = "always")`.
  freq <- match.arg(
    getOption("vista.deprecation_frequency", "regularly"),
    c("regularly", "once", "always")
  )
  if (identical(freq, "always")) {
    rlang::warn(message = msg, class = cls)
  } else {
    rlang::warn(
      message = msg,
      class = cls,
      .frequency = freq,
      .frequency_id = paste(fun, old, sep = "/")
    )
  }

  transform(value)
}

#' Signal a defunct VISTA argument
#'
#' Companion to [.vista_deprecate_arg()] used once an alias graduates from
#' warning to error. Not wired up yet; see `?VISTA-deprecated` for timelines.
#'
#' @inheritParams .vista_deprecate_arg
#' @return Never returns; always aborts.
#' @keywords internal
#' @noRd
.vista_defunct_arg <- function(old, new = NULL, fun) {
  row <- .vista_deprecation_row(fun, old)
  if (is.null(new) && !is.null(row)) new <- row$new_arg
  target <- if (!is.null(new) && nzchar(new)) {
    sprintf("Use `%s` instead.", new)
  } else {
    "It no longer has any effect and must be removed from your call."
  }
  rlang::abort(
    message = sprintf("`%s` in `%s()` is defunct. %s", old, fun, target),
    class = c("vista_defunct_arg", "defunctError")
  )
}

#' @keywords internal
#' @noRd
.vista_did_you_mean <- function(bad, candidates, n = 1L) {
  if (!length(candidates)) return(character(0))
  d <- utils::adist(bad, candidates, ignore.case = TRUE)[1, ]
  ok <- d <= pmax(2L, floor(nchar(bad) / 3))
  if (!any(ok)) return(character(0))
  utils::head(candidates[order(d)][ok[order(d)]], n)
}

#' Reject unknown arguments passed through `...`
#'
#' Generalizes the guard previously hand-written in [get_pathway_heatmap()].
#' VISTA plotting functions forward `...` to third-party plotting engines, which
#' silently ignore names they do not recognise — so a typo such as
#' `gene = my_genes` would previously fall through and plot the default gene set
#' instead of erroring.
#'
#' @param dots A list, typically `list(...)`.
#' @param fun Character scalar naming the calling function.
#' @param allowed Character vector of argument names that are legitimately
#'   forwarded (for example the formals of the downstream plotting function).
#' @param blocked Character vector of names the caller manages itself and which
#'   must never be forwarded.
#'
#' @return `invisible(TRUE)`; called for its side effect.
#' @keywords internal
#' @noRd
.vista_check_dots <- function(dots, fun, allowed = character(0), blocked = character(0)) {
  nms <- names(dots)
  if (is.null(nms) || !length(nms)) return(invisible(TRUE))

  unnamed <- !nzchar(nms)
  if (any(unnamed)) {
    cli::cli_abort(
      "{.fun {fun}} received {sum(unnamed)} unnamed argument{?s} in {.arg ...}. All {.arg ...} arguments must be named."
    )
  }

  hit_blocked <- intersect(nms, blocked)
  if (length(hit_blocked)) {
    cli::cli_abort(
      "{cli::qty(hit_blocked)}Argument{?s} {.val {hit_blocked}} {?is/are} managed by {.fun {fun}} and cannot be passed through {.arg ...}."
    )
  }

  known <- union(allowed, names(formals(fun, envir = asNamespace("VISTA"))))
  unknown <- setdiff(nms, known)
  if (!length(unknown)) return(invisible(TRUE))

  bullets <- vapply(unknown, function(nm) {
    sug <- .vista_did_you_mean(nm, known)
    if (length(sug)) sprintf("%s (did you mean `%s`?)", nm, sug[[1]]) else nm
  }, character(1))

  cli::cli_abort(c(
    "{cli::qty(unknown)}{.fun {fun}} received unknown argument{?s} in {.arg ...}.",
    stats::setNames(bullets, rep("x", length(bullets)))
  ))
}

#' Deprecated and defunct arguments in VISTA
#'
#' @description
#' VISTA follows the Bioconductor deprecation contract: a renamed argument warns
#' for one release, becomes defunct in the next, and is removed in the one after.
#' Deprecated arguments keep working exactly as before while they warn.
#'
#' @details
#' The table below is generated from the package's internal registry, so it is
#' always in sync with the code. Warnings are classed
#' `c("vista_deprecated_arg", "deprecatedWarning")` and are emitted at most once
#' per session per argument, so they will not flood a loop.
#'
#' To silence them while you migrate:
#' \preformatted{
#' suppressWarnings(get_corr_heatmap(v, show_corr_values = TRUE))
#' }
#'
#' To find every deprecated call in your own scripts, promote the warnings to
#' errors:
#' \preformatted{
#' options(warn = 2)
#' }
#'
#' @section Registry:
#' \Sexpr[stage=render,results=rd]{VISTA:::.vista_deprecation_rd()}
#'
#' @name VISTA-deprecated
#' @aliases VISTA-deprecated
#' @return Not applicable; this page documents argument deprecations only.
#' @examples
#' # Inspect the timeline for a specific function
#' reg <- VISTA:::.vista_deprecations()
#' reg[reg$fun == "get_corr_heatmap", c("old_arg", "new_arg", "defunct_in")]
NULL

#' @keywords internal
#' @noRd
.vista_deprecation_rd <- function() {
  reg <- .vista_deprecations()
  reg <- reg[order(reg$fun, reg$old_arg), , drop = FALSE]
  esc <- function(x) gsub("([\\\\{}%])", "\\\\\\1", x)
  body <- paste0(
    "\\code{", esc(reg$fun), "()} \\tab \\code{", esc(reg$old_arg), "} \\tab ",
    ifelse(nzchar(reg$new_arg), paste0("\\code{", esc(reg$new_arg), "}"), "removed"),
    " \\tab ", esc(reg$defunct_in), " \\cr",
    collapse = "\n"
  )
  paste0(
    "\\tabular{llll}{\n",
    "\\strong{Function} \\tab \\strong{Deprecated} \\tab \\strong{Use instead} \\tab \\strong{Defunct in} \\cr\n",
    body, "\n}"
  )
}
