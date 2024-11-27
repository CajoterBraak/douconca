#' @noRd
#' @keywords internal
f_inertia <- function(object) {
  # function f_inertia uses vegan 2.6-4 internal structure
  # object: a dccav object, results of dc_CA
  # value:  a matrix (currently with 1 column) with in first column the inertias
  # if (is.null(object$data$Y)) {
  #   object$env_explain <- object$inertia["env_explain"]
  # } else {
  #   ms <- f_wmean(object$formulaTraits, tY = object$data$Y, object$data$dataTraits,
  #                 weights=object$weights, name= "CWM")
  #   object$env_explain <- ms$explained
  # }
  if (is.null(object$CCAonTraits)) {
    total <- object$tot.chi
    conditionT <- object$conditionT
    env_explain <- object$inertia["env_explain"]
  } else {
    total <- object$CCAonTraits$tot.chi
    conditionT <- object$CCAonTraits$pCCA$tot.chi
    mt <- try(f_wmean(object$formulaEnv, tY = t(object$data$Y)/sum(object$data$Y), object$data$dataEnv,
                  weights=object$weights, name= "SNC"))
    if (inherits(mt, "try-error")) {
      warning("singular environment data. Env_explain not available.\n") 
      env_explain <- NA
    } else {env_explain <- mt$explained; names(env_explain) <- NULL}
  }
  if (is.na(env_explain)) env_explain <- NULL
  inertia <- cbind(c(total = total,
                     conditionT = conditionT,
                     traits_explain = object$RDAonEnv$tot.chi,
                     env_explain = env_explain,
                     conditionTE = object$RDAonEnv$pCCA$tot.chi,
                     constraintsTE = object$RDAonEnv$CCA$tot.chi))
  colnames(inertia) <- "weighted variance"
  expla <- c("total inertia (= weighted variation)",
             "variation fitted by the trait condition", 
             "trait-constrained variation", 
             "environment-constrained variation",
             "trait-constrained variation explained by the condition in formulaEnv",
             "trait-constrained variation explained by the predictors in formulaEnv")
  names(expla) <- c("total", "conditionT", "traits_explain", "env_explain",
                    "conditionTE", "constraintsTE")
  attr(inertia, which = "meaning") <- 
    matrix(expla[rownames(inertia)], ncol = 1, 
           dimnames = list(rownames(inertia), "meaning"))
  return(inertia)
}

#' @noRd
#' @keywords internal
get_QR <- function(object, model = "CCA"){
  # function get_QR uses vegan 2.6-4 internal structure
  # gets the qr decompostion of object
  # model = "CCA" or "pCCA"
  if (model == "CCA") {
    QR <- object$CCA$QR 
  } else if (model == "pCCA") {
    QR <- object$pCCA$QR
  } else {
    stop("model not supported.\n")
  }
  return(QR)
}
