#' Print a summary of a wrda object
#' 
#' @param x a wrda object from \code{\link{wrda}}
#' @param ...  Other arguments passed to the function (currently ignored).
#' 
#' @return No return value, results are printed to console.
#' 
#' @example demo/dune_wrda.r
#'
#' @export
print.wrda <- function(x, 
                       ...) {
  print_dcca(x)
}
