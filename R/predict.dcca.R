#' @title Prediction for double-constrained correspondence analysis (dc-CA)
#'
#' @description
#' Prediction of traits from environment, environment from traits and response 
#' from trait and environment data.
#' 
#' With \code{type = "traits"} and \code{newdata = NULL}, predict gives the 
#' fitted mean traits, \emph{i.e.} the fitted community weighted means.
#' With \code{type = "env"} and \code{newdata = NULL}, predict gives the 
#' fitted mean environment, \emph{i.e.} the fitted species niche centroids.
#'
#' @param object return value of \code{\link{dc_CA}}.
#' @param ...  Other arguments passed to the function (currently ignored).
#' @param type type of prediction, \code{c("env", "traits", "response", 
#' "reg_env", "reg_trait")} for environmental values, values of traits, 
#' response (expected abundance), regression coefficients for environmental 
#' and trait predictors.
#' @param newdata Data in which to look for variables with which to predict.
#' For \code{type = "reg_env" or "reg_trait"} newdata is ignored. For 
#' \code{type = "env" or "trait"}, \code{newdata} is a data frame of trait and 
#' environmental values, respectively, which are used for the prediction. If 
#' omitted, fitted values are generated.
#' For \code{type = "response"}, newdata is a list of two data frames with 
#' trait and environmental values in this order, \emph{e.g.} 
#' \code{list(traits = dataTraits, env = dataEnv)}.
#' @param rank rank or number of axes to use. Default "full" for all axes 
#' (no rank-reduction).
#' 
#' @details
#' Variables that are in the model but not in \code{newdata} are set to their 
#' weighted means in the training data. Predictions are thus at the (weighted)
#' mean of the quantitative variables not included. Predictions with 
#' not-included factors are at the weighted mean (none of the factor effects 
#' are included).
#'
#' For \code{type = "response"} and non-null newdata, the species weights of 
#' the training are used; the site weights are taken equal. Many of the 
#' predicted values may be negative, indicating expected absences (0) or small
#' expected response values.
#' 
#' Regression coefficients obtained with \code{type = "reg_env"} or  
#' \code{type = "reg_traits"} are for standardized traits and environmental
#' variables.
#' 
#' @return a matrix with the predictions. The exact content of the matrix 
#' depends on the \code{type} of predictions that are being made.
#'
#' @example demo/dune_dcCA_predict.R
#' 
#' @export
predict.dcca <- function(object,
                         ...,
                         type = c("env", "traits", "response", 
                                  "reg_env", "reg_traits"),
                         rank = "full",
                         newdata = NULL) {
  type <- match.arg(type)
  if (rank == "full") {
    rank <- length(object$eigenvalues)
  }
  if (type == "response") {
    if (is.null(newdata)) {
      newdata <- list(NULL, NULL)
    }
    newdata1 <- list(
      # env prediction requires trait data
      traits = check_newdata(object, newdata[[1]], "env"), 
      # trait prediction requires env data
      env = check_newdata(object, newdata[[2]], "traits") 
    )
  } else if (type %in% c("env", "traits")) {
    newdata1 <- check_newdata(object, newdata, type)
  }
  ret <- switch(type,
                env = predict_env(object, newdata1, rank),
                traits = predict_traits(object, newdata1, rank),
                response = predict_response(object, newdata1, rank),
                reg_env = predict_regr_env(object, rank),
                reg_traits = predict_regr_traits(object, rank)
                
  )
  return(ret)
}
