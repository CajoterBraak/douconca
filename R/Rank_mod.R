#' @param object result of dc_CA (class dcca or dccav)
#' @param partial logical (TRUE not implemented)
#' 
#' @noRd
#' @keywords internal
Rank_mod <- function(object, 
                     partial = FALSE) {
  # returns the rank of a dc-CA model
  if (inherits(object, c("dccav", "dcca", "wrda"))) {
    rr <- length(object$eigenvalues)
  } else if (inherits(object, "cca") && !partial) {
    rr <- length(vegan::eigenvals(object, model = "constrained"))
  } else if (inherits(object, "cca") && partial) {
    stop ("rank of pCCA model not implemented")
  } else {
    stop("object in Rank_mod must be of class cca or dccav")
  }
  return(rr)
}
