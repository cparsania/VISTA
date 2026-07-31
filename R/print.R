#' Print a VISTA object
#'
#' S3 entry point kept for code that calls `print()` explicitly (including
#' knitr, which prints rather than auto-displaying). It forwards to the S4
#' `show()` method so both routes produce identical output.
#'
#' @param x A VISTA object.
#' @param ... Ignored.
#'
#' @return `x`, invisibly. Called for its output.
#' @examples
#' v <- example_vista()
#' print(v)
#' @export
#' @method print VISTA
print.VISTA <- function(x, ...) {
  methods::show(x)
  invisible(x)
}
