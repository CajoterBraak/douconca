#' @title Fitted values of double-constrained correspondence analysis (dc-CA)
#'
#' @description
#' Community weighted means (CWM) and species-niche centroids (SNC),
#' as fitted (in full or reduced rank) from the environmental data and
#' trait data, respectively, and  
#' the fitted response from trait and environment data.
#'
#' @param object return value of \code{\link{dc_CA}}.
#' @param ...  Other arguments passed to the function (currently ignored).
#' @param type type of prediction, \code{c( "CWM","SNC", "response")} 
#' for environmental values, values of traits, 
#' response (expected abundance).
#' @param rank rank (number of axes to use). Default "full" for all axes 
#' (no rank-reduction).
#' 
#' @details
#'
#' If \code{type="response"} the rowsums of \code{object$data$Y} are used
#' to scale the fit to these sums.
#' Many of the predicted response values may be negative, 
#' indicating expected absences (0) or small expected response values.
#' 
#' @return a matrix with fitted value. The exact content of the matrix 
#' depends on the \code{type} of fits that are asked for.
#' 
#'
#' @example demo/dune_dcCA_fitted.R
#' 
#' @export
fitted.dcca <- function(object,
                         ...,
                         type = c("CWM","SNC", "response"),
                         rank = "full") {
  newdata <- NULL
  type <- match.arg(type)
  if (rank == "full") rank <- length(object$eigenvalues)
  if (type == "response") {
      newdata1 <- list(
      # env prediction requires trait data
      traits = object$data$dataTraits, 
      # trait prediction requires env data
      env = object$data$dataEnv) 
  } 

  ret <- switch(type,
                SNC      = predict_env(object,  object$data$dataTraits, rank),
                CWM      = predict_traits(object, object$data$dataEnv, rank),
                response = predict_response(object, newdata1, rank, object$weights)
  )
  if (type =="response"){
    if (is.null(object$data$Y)) totsum <- 1 else totsum <- sum(object$data$Y) 
    ret <- ret * totsum
  }
  return(ret)
}
