#'
#' @param object result of dc_CA_vegan of class dccav
#' @param partial logical (TRUE not implemented)
#' @noRd
# @export
Rank_mod <- function(object, partial = FALSE){
  # returns the rank of a dc-CA model
  if ("dccav" %in% class(object)) {
    rr <- length(object$eigenvalues)
  } else if ("cca" %in% class(object) && !partial){
    rr <- length(vegan::eigenvals(object, model = "constrained"))
  } else if ("cca" %in% class(object) && partial){
    stop ("rank of pCCA model not implemented")
    #rr <- get_QR(object, model = "pCCA")$rank
  } else stop("object in Rank_mod must be of class cca or dccav")
  return(rr)
}
