#' @title Coefficients of double-constrained correspondence analysis (dc-CA)
#'
#' @description
#' Fourth-corner coefficients and regression coefficients 
#' (of full or reduced rank) to predict traits from environment, 
#' environment from traits and response from trait and environment data.
#'
#' @param object return value of \code{\link{dc_CA}}.
#' @param ...  Other arguments passed to the function (currently ignored).
#' @param type type of coefficients, 
#' \code{c("fourth_corner", "all_reg", "env2traits_reg", "traits2env_reg")} for
#' fourth-corner coefficients and regression coefficients for
#' all trait x environmental predictors, environmental predictors only 
#' and trait predictors only for prediction of the (transformed)
#' response, traits and environmental values, respectively.
#' 
#' @details
#' 
#' Regression coefficients are for standardized traits and environmental variables.
#' 
#' With covariates, \code{coef()} gives partialfourth-corner correlations. 
#' With \code{rank = 2}, \code{coef()} gives the two-dimensional approximation
#' of the full-rank fourth-corner correlations in the biplot that displays the
#' traits and environmental variables at arrow heads or points 
#' at \code{scores(mod, display = c("bp", "bp_traits"))}.
#' 
#' 
#' @return a matrix with coefficients. The exact content of the matrix 
#' depends on the \code{type} of coefficient that are asked for.
#' 
#' Regression coefficients for a response variable 
#' are usually column-vectors. 
#' With \strong{X} the matrix of units-by-predictors
#' and \strong{B} the matrix of predictors-by-response-variables, 
#' predictions or fits are of the form \strong{Y = XB}.
#' Analogously, \code{type = "trait2env"} gives a trait-by-environment matrix and
#' \code{type = "env2traits"} gives an environment-by-trait matrix.
#' 
#'
#' @example demo/dune_dcCA_coef.R
#' 
#' @export
coef.dcca <- function(object,
                         ...,
                         type = c("fourth_corner", "all_reg", "env2traits_reg", "traits2env_reg"),
                         rank = "full",
                         newdata = NULL) {
  type <- match.arg(type)
  if (rank == "full") {
    rank <- length(object$eigenvalues)
  }
  ret <- switch(type,
                fourth_corner  = predict_fc(object, rank),
                all_reg        = predict_regr_all(object, rank),
                env2traits_reg = predict_regr_env(object, rank),
                traits2env_reg = predict_regr_traits(object, rank)
  )
  return(ret)
}
