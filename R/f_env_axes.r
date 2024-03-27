#' Utility function: obtain column (sites) related scores
#' @param out object from \code{\link{dc_CA_vegan}}
#' @param which_cor character or names of environmental variables
#' for which inter-set correlations must calculated.
#' Default: "in_model" for all environmental variables in the model specified by \code{formulaEnv}.
# @noRd
# @export
f_env_axes <- function(out, which_cor = "in model"){
  # which_cor character or names of environmental variables
  # for which inter-set correlations must calculated.
  #  Default "in_model" for all environmental variables in the model

  #print(names(out))
  myconst <- sqrt(out$Nobs*out$RDAonEnv$tot.chi)
  CWM <- as.matrix(vegan::scores(out$RDAonEnv, display = c("sites"), scaling =  "sites",
                choices =  seq_len(Rank_mod(out$RDAonEnv)), const = myconst))
  res <- calculate_b_se_tval(get_QR(out$RDAonEnv), y=CWM,
                             w = out$weights$rows,  scale2 = 0, name = "CWM")
  c_env_normed <- res$coef_normed
  attr(c_env_normed, which = "warning") <-"The t-values are optimistic, i.e. an underestimate of their true absolute value"


  # correlations of the dataEnv with the CWMs wrt the  axes
  if (which_cor == "in model"){
    #in_model <- colnames(out$data$dataEnv)%in% rownames(attr(stats::terms(out$RDAonEnv), which = "factors"))
     in_model <- get_focal_and_conditioning_factors(out$RDAonEnv, factors_only = FALSE)$`focal factor`
    } else in_model = which_cor
  env0 <-  stats::model.matrix(~.-1, constrasts = FALSE, data = out$data$dataEnv[, in_model, drop= FALSE])
  Cor_Env_CWM <- wcor(env0, CWM, w = out$weights$rows)
  colnames(Cor_Env_CWM) <- paste("CWM-ax", seq_len(ncol(Cor_Env_CWM)), sep= "")
  attr(Cor_Env_CWM, which = "meaning")<- "inter set correlation, correlation between environmental variables and CWM of axes"


  out2 <- list(site_scores = list(site_scores_unconstrained = res$y,
                                     lc_env_scores = res$fitted), c_env_normed= c_env_normed, b_se= res$b_se, R2_env = res$R2, correlation = Cor_Env_CWM)
  return(out2)
}
