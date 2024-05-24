#' @title Prediction for double-constrained correspondence analysis (dc-CA)
#'
#' @description
#' Prediction of traits from environment and environment from traits
#'
#' With \code{type = "traits"} and \code{newdata=NULL}, predict gives the fitted mean traits,
#' \emph{i.e.} the fitted community weighted means.
#' With \code{type = "env"} and \code{newdata=NULL}, predict gives the fitted mean environment,
#' \emph{i.e.} the fitted species niche centroids.
#'
#' @param object return value of \code{\link{dc_CA}}.
#' @param type type of prediction, c("env", "traits","response","reg_env","reg_trait") for
#' environmental values, values of traits,
#' regression coefficients for environmenal and trait predictors.
#' @param newdata An obtional data frame in which to look for variable with which to predict. If omitted,
#' the fitted values are generated.
#' @param rank rank or number of axes to use. Default "full" for all axes (no rank-reduction).
#' @details
#' Variables that are in the model but not in \code{newdata} are set to their weighted means in
#' the training data.
#' Predictions are thus at the (weighted) mean of the quantitative variables not included.
#' Predictions with not-included factors are at the weighted mean (none of the factor effects are included).
#'
#' @example demo/dune_dcCA_predict.R
#' @export

predict.dcca <- function(object,
                         type = c("env", "traits","response","reg_env","reg_traits"),
                         rank = "full",
                         newdata=NULL){
  type <- match.arg(type)
  if (rank == "full") rank = length(object$eigenvalues)
  if (type %in% c("traits","regr_env")) {
    ff <- get_Z_X_XZ_formula(object$formulaEnv)
    c_normed <- object$c_env_normed
    if (is.null(newdata)) newdata <- object$data$dataEnv
  } else {
    ff <- get_Z_X_XZ_formula(object$formulaTraits)
    c_normed <- object$c_traits_normed
    if (is.null(newdata)) newdata <- object$data$dataTraits
  }

  if (type %in% c("traits","env")) {
      nams <- !ff$all_nams %in% names(newdata)
      nams <-ff$all_nams[nams]
      for (n in nams){
          if (n %in% rownames(c_normed)) newdata[[n]] <- c_normed[n,"Avg"] # set to average
          else newdata[[n]] <- 0 # for factors
      }
      newdata1 <- stats::model.matrix(ff$formula_XZ, constrasts = FALSE, data = newdata)[,-1,drop=FALSE]

  }
  ret <- switch(type,
                env = predict_env(object, newdata1, rank),
                traits = predict_traits(object,newdata=newdata1, rank),
                reg_env = predict_regr_env(object, rank),
                reg_traits = predict_regr_traits(object, rank)

                )
 return(ret)
}
