#' Utility function: obtain row (sites) related scores
#' @param out object from \code{\link{dc_CA}}, \code{\link[vegan]{cca}} or  \code{\link[vegan]{rda}}
#' @param which_cor character or names of environmental variables
#' @param wrda logical. When FALSE out is an object from \code{\link[vegan]{rda}} or \code{\link[vegan]{cca}}, else from \code{\link{wrda}}
#' for which inter-set correlations must calculated.
#' Default: "in_model" for all environmental variables in the model specified by \code{formulaEnv}.
#' @noRd
# @export
f_env_axes <- function(out, which_cor = "in model"){
  # which_cor character or names of environmental variables
  # for which inter-set correlations must calculated.
  #  Default "in_model" for all environmental variables in the model

  #print(names(out))
  if ("dcca" == class(out)[1]) {
  dataEnv <- out$data$dataEnv
  formulaEnv <- out$formulaEnv
  out <- out$RDAonEnv
  }

  w <- out$weights$rows

  if ("dccav" %in% class(out)) {

    myconst <- sqrt(out$Nobs*out$RDAonEnv$tot.chi)
    CWM <- as.matrix(vegan::scores(out$RDAonEnv, display = c("sites"), scaling =  "sites",
                          choices =  seq_len(Rank_mod(out$RDAonEnv)), const = myconst))
    dataEnv <- out$data$dataEnv;
    formulaEnv <-out$formulaEnv
    out <- out$RDAonEnv;
  } else if ("wrda" %in% class(out)) {
    CWM <- out$site_axes$site_scores$site_scores_unconstrained
    dataEnv <- out$data
    formulaEnv <- out$formula
  } else stop("input must be a result of dc_CA or wrda")
  QR <- get_QR(out)
  res <- calculate_b_se_tval(QR, y=CWM,
                             w = w,  scale2 = 0, name = "CWM")
  c_env_normed <- res$coef_normed
  attr(c_env_normed, which = "warning") <-"The t-values are optimistic, i.e. an underestimate of their true absolute value"


  # correlations of the dataEnv with the CWMs wrt the  axes
  if (which_cor[1] == "in model"){
      fX <- get_Z_X_XZ_formula(formulaEnv)$formula_X0
      env0 <-  stats::model.matrix(fX, constrasts = FALSE, data = dataEnv)
  } else {
    whichc = which_cor
    env0 <-  stats::model.matrix(~.-1, constrasts = FALSE, data = dataEnv[, whichc, drop= FALSE])
    }
  Cor_Env_CWM <- wcor(env0, CWM, w = w)
  colnames(Cor_Env_CWM) <- paste("CWM-ax", seq_len(ncol(Cor_Env_CWM)), sep= "")
  attr(Cor_Env_CWM, which = "meaning")<- "inter set correlation, correlation between environmental variables and CWM of axes"


  out2 <- list(site_scores = list(site_scores_unconstrained = res$y,
                                     lc_env_scores = res$fitted), c_env_normed= c_env_normed, b_se= res$b_se, R2_env = res$R2, correlation = Cor_Env_CWM)
  return(out2)
}
