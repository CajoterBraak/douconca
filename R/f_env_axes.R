#' Utility function: obtain row (sites) related scores
#' 
#' @param out object from \code{\link{dc_CA}}, \code{\link[vegan]{cca}} or  
#' \code{\link[vegan]{rda}}
#' @param which_cor character or names of environmental variables
#' @param wrda logical. When FALSE out is an object from 
#' \code{\link[vegan]{rda}} or \code{\link[vegan]{cca}}, else from 
#' \code{\link{wrda}} for which inter-set correlations must calculated.
#' Default: "in_model" for all environmental variables in the model specified
#' by \code{formulaEnv}.
#' 
#' @noRd
#' @keywords internal
f_env_axes <- function(out, 
                       which_cor = "in model") {
  # which_cor character or names of environmental variables
  # for which inter-set correlations must calculated.
  #  Default "in_model" for all environmental variables in the model
  if (inherits(out, "dcca", which = TRUE) == 1) {
    dataEnv <- out$data$dataEnv
    formulaEnv <- out$formulaEnv
    out <- out$RDAonEnv
  }
  w <- out$weights$rows
  if (inherits(out, "dccav")) {
    myconst <- sqrt(out$Nobs * out$RDAonEnv$tot.chi)
    CWM <- as.matrix(scores(out$RDAonEnv, display = "sites", 
                            scaling = "sites",
                            choices = seq_len(rank_mod(out$RDAonEnv)), 
                            const = myconst))
    dataEnv <- out$data$dataEnv
    formulaEnv <-out$formulaEnv
    out <- out$RDAonEnv
  } else if (inherits(out, "wrda")) {
    CWM <- out$site_axes$site_scores$site_scores_unconstrained
    dataEnv <- out$data
    formulaEnv <- out$formula
  } else {
    stop("input must be a result of dc_CA or wrda.\n")
  }
  formulaEnv <- change_reponse(formulaEnv, "Y", dataEnv)
  QR <- get_QR(out)
  res <- calculate_b_se_tval(QR, y = CWM, w = w, scale2 = 0, name = "CWM")
  c_env_normed <- res$coef_normed
  # NB: if vegan::rda, then SDS is a factor sqrt(Nobs) too large
  if (diff(range(w)) < 1.0e-4 / length(w)) 
    c_env_normed[, 2] <-c_env_normed[, 2]/sqrt(length(w))
  attr(c_env_normed, which = "warning") <-
    "The t-values are optimistic, i.e. an underestimate of their true absolute value"
  # NB: the
  #fX <- get_Z_X_XZ_formula(formulaEnv, dataEnv)$formula_XZ
  #envXZ <- model.matrix(formula = fX, data = dataEnv)
  #msdXZ <- mean_sd_w(envXZ, w = w)
  #if (abs(c_env_normed[, 2] - msdXZ$sd[1,]>0.0006)) stop("SDS wrong")
  ###NB end

      
  
  # correlations of the dataEnv with the CWMs wrt the  axes
  if (which_cor[1] == "in model"){
    fX <- get_Z_X_XZ_formula(formulaEnv, dataEnv)$formula_X1
  } else {
    whichc <- which_cor
    fX <- as.formula(paste("~", paste0(whichc, collapse = "+")))
  }
  env0 <- modelmatrixI(formula = fX, data = dataEnv, XZ = FALSE)
  msd <- mean_sd_w(env0, w = w)
  #nams <- intersect(rownames(c_env_normed), colnames(msd$sd))
  #c_env_normed[nams,"SDS"] <- msd$sd[1, nams]
  Cor_Env_CWM <- wcor(env0, CWM, w = w)
  colnames(Cor_Env_CWM) <- paste0("CWM-ax", seq_len(ncol(Cor_Env_CWM)))
  attr(Cor_Env_CWM, which = "meaning") <- 
    "inter set correlation, correlation between environmental variables and CWM of axes"
  out2 <- list(site_scores = list(site_scores_unconstrained = res$y,
                                  lc_env_scores = res$fitted), 
               c_env_normed = c_env_normed, 
               b_se = res$b_se, msd = msd,
               R2_env = res$R2, 
               correlation = Cor_Env_CWM)
  return(out2)
}
