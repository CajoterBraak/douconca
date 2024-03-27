#' Takes the scores function to the dcCA library
#'
#' @param x in package \code{dcCA}, an object of \code{\link{dc_CA_vegan}}, in package \code{vegan}
#' an  ordination object, see also \code{\link[vegan]{cca.object}}.
#' @param ...  Other arguments passed to the function.
#' @seealso \code{\link{scores.dccav}} and \code{\link[vegan]{scores}}
#' @export
"scores" <-
function(x, ...) UseMethod("scores")

