#' Print a summary of a dc-CA object.
#' 
#' @param x a dc-CA object from \code{\link{dc_CA}}.
#' @param ...  Other arguments passed to the function (currently ignored).
#' 
#' @details
#' \code{x <- print(x)} is more efficient for \code{\link{scores.dcca}} than 
#' just \code{print(x)} if \code{\link{dc_CA}} is called with 
#' \code{verbose = FALSE}).
#' 
#' @return No return value, results are printed to console.
#'
#' @example demo/dune_dcCA.r
#'
#' @export
print.dcca <- function(x,
                       ...) {
  print_dcca(x)
}
