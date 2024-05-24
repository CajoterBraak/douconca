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
  } else {
    ff <- get_Z_X_XZ_formula(object$formulaTraits)
  }

  if (type %in% c("traits","env")) {
    if (is.null(newdata)){
      if (type =="traits") newdata <- object$data$dataEnv else #if (type == "env")
        newdata <- object$data$dataTraits
    }
      nams <- !ff$all_nams %in% names(newdata)
      nams <-ff$all_nams[nams]
      for (n in nams){
        # ff.data <- get_Z_X_XZ_formula(object$formulaEnv,object$data$dataEnv)
        # if (n %in% c(ff.data$focal_factor, ff.data$Condi_factor)){
        #   #n is a factor: todo, this is a hack
        #   newdata[[n]] <- 0
        # } else
          newdata[[n]] <- 0
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
