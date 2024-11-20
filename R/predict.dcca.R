#' @title Prediction for double-constrained correspondence analysis (dc-CA)
#'
#' @description
#' Prediction of traits from environment, environment from traits and response 
#' from trait and environment data.
#'
#' @param object return value of \code{\link{dc_CA}}.
#' @param ...  Other arguments passed to the function (currently ignored).
#' @param type type of prediction, \code{c("envFromTraits", "traitsFromEnv",
#' "response", "lc", "lc_traits")} 
#' for environmental values, values of traits, 
#' response (expected abundance) and constrained scores for sites and species. 
#' \code{"SNC"} is equivalent with \code{c("envFromTraits")}.
#' \code{"CWM"} is equivalent with \code{c("traitsFromEnv")}.

#' @param newdata Data in which to look for variables with which to predict.
#' For \code{type = "envFromTraits" or "traitsFromEnv"} or 
#'  \code{type = "lc_traits" or "lc"},
#' \code{newdata} is a data frame of trait and environmental values, respectively, 
#' which are used for the prediction or the calculation of scores. 
#' If omitted, fitted values are generated (use \code{\link{fitted.dcca}} instead).
#' For \code{type = "response"}, newdata is a list of two data frames with 
#' trait and environmental values in this order, \emph{e.g.} 
#' \code{list(traits = dataTraits, env = dataEnv)}.
#' @param rank rank (number of axes to use). Default "full" for all axes 
#' (no rank-reduction).
#' @param weights list of weights of species and of sites in \code{newdata} when
#' \code{type = "response"}, else ignored (default NULL
#' yielding equal species and site weights, both summing to 1). 
#' Example: weights = list(species = c(100,1,1), sites = c(1,1,1,1)), in that
#' order, with traits of three new species in newdata[[1]] and 
#' environmental values (and levels of factors) of four new sites in newdata[[2]]. 
#' Species weights are scaled to a sum of one.
#' @inheritParams scores.dcca
#' 
#' @details
#' Variables that are in the model but not in \code{newdata} are set to their 
#' weighted means in the training data. Predictions are thus at the (weighted)
#' mean of the quantitative variables not included. Predictions with 
#' not-included factors are at the reference level (the first level of the factor).
#'
#' For \code{type = "response"} and non-null \code{newdata}, the species weights of 
#' the training are used; the site weights are taken equal. Many of the 
#' predicted values may be negative, indicating expected absences (0) or small
#' expected response values.
#' 
#' With \code{type = "traitsFromEnv"} and \code{newdata = NULL}, predict gives the 
#' fitted mean traits, \emph{i.e.} the fitted community weighted means.
#' With \code{type = "envFromTraits"} and \code{newdata = NULL}, predict gives the 
#' fitted mean environment, \emph{i.e.} the fitted species niche centroids
#' (see \code{\link{fitted.dcca}}). See \code{fitted.dcca}.
#' 
#' @return a matrix with the predictions. The exact content of the matrix 
#' depends on the \code{type} of predictions that are being made.
#'
#' @example demo/dune_dcCA_predict.R
#' 
#' @export
predict.dcca <- function(object,
                         ...,
                         type = c("envFromTraits", "traitsFromEnv", "response", 
                                  "SNC", "CWM","lc","lc_traits"),
                         rank = "full",
                         newdata = NULL, weights = NULL,
                         scaling = "symmetric"
                         ) {
  type <- match.arg(type)
  if (rank == "full") rank <- length(object$eigenvalues)
  if (type == "response") {
    if (is.null(newdata)) {
      newdata <- list(traits = object$data$dataTraits, env = object$data$dataEnv)
      weights <- object$weights
    } else {
      if (is.list(newdata)&& length(newdata)==2) {
        newdata <- list(
      # env prediction requires trait data
        traits = newdata[[1]], 
      # trait prediction requires env data
        env =  newdata[[2]])
      } else stop("For prediction of response," ,
              " newdata must be a list of trait and env data")
      if (is.null(weights)) weights = list(species = NULL, sites = NULL)
      if (is.null(weights[[1]])) weights$species <- 
          rep(1/nrow(newdata[["traits"]]), nrow(newdata[["traits"]]))
      if (is.null(weights[[2]])) weights$sites  <- 
          rep(1/nrow(newdata[["env"]]),nrow(newdata[["env"]]))
      if (!length(weights[[1]])== nrow(newdata[["traits"]])){
        weights[[1]] <- 
          rep(1/nrow(newdata[["traits"]]), nrow(newdata[["traits"]]))
        warning("length of weights for species does not match new trait data. ",
                "Species weights reset to equal weights.\n")
      }
      if (!length(weights[[2]])== nrow(newdata[["env"]])){
        weights[[2]] <- 
          rep(1/nrow(newdata[["env"]]), nrow(newdata[["env"]]))
        warning("length of weights for sites does not match new environment data. ",
                "Site weights reset to equal weights.\n")
      }
    }
  } else {
    if (is.null(newdata)) {
      if (type %in% c("traitsFromEnv","SNC") )newdata <- object$data$dataEnv else
        newdata <- object$data$dataTraits
    }
  }

  ret <- switch(type,
                envFromTraits = predict_env(object, newdata, rank),
                traitsFromEnv = predict_traits(object, newdata, rank),
                SNC = predict_env(object, newdata, rank),
                CWM = predict_traits(object, newdata, rank),
                response = predict_response(object, newdata, rank, weights),
                lc = predict_lc(object, newdata, rank, scaling = scaling),
                lc_traits = predict_lc_traits(object, newdata, rank, 
                                              scaling = scaling)
  )
  return(ret)
}
