#' Print a VISTA object like a SummarizedExperiment
#'
#' Forwards to SummarizedExperiment's `show()` so the output is identical to a
#' plain SE. Invisibly returns `x`.
#'
#' @param x A VISTA object.
#' @param ... Ignored.
#'
#' @export
#' @method print VISTA
print.VISTA <- function(x, ...) {
  se_show <- methods::selectMethod("show", "SummarizedExperiment")
  se_show(as(x, "SummarizedExperiment"))
  invisible(x)
}

#' @export
#' @method print vista
print.vista <- function(x, ...) {
  # Back-compat if anything in your code still calls print() on old 'vista' class
  print.VISTA(x, ...)
}
