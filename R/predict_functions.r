predict_regr_env <- function(object, rank){
 sc <- scores(object, choices = seq_len(rank), display = c("reg", "bp_traits"))
 B_env_regr <-sc$regression
 C_trait_bip <- sc$biplot_traits
 reg <-   B_env_regr[,-c(1,2,3), drop = FALSE]%*%  t(C_trait_bip)
 return(reg)
}

predict_regr_traits <- function(object, rank){
  sc <- scores(object, choices = seq_len(rank), display = c("reg_traits", "bp"))
  B_traits_regr <-sc$regression_traits
  C_env_bip <- sc$biplot
  reg <-   B_traits_regr[,-c(1,2,3), drop = FALSE]%*%  t(C_env_bip)
  return(reg)
}

predict_traits <- function(object,  newdata, rank){
  # missing factors in newdata : no contribution of the missing factors
  reg <- predict_regr_env(object, rank)
  nms <- rownames(reg)
  nms <- nms[nms %in% colnames(newdata)] # still needed for missing factors; see predict.dcca.r
  dat0 <- newdata[,nms, drop = FALSE]
  dat <- scale_data(dat0, mean_sd = object$c_env_normed[,c(1,2)])
  #pred_scaled <- dat %*% reg[colnames(dat),, drop = FALSE ]
  # reg may not have the covariates (e.g. when dc-CA was started from CWM) and
  # newdata may not have all predictors
  nams <- intersect(colnames(dat), rownames(reg))
  pred_scaled <- dat[, nams,drop = FALSE] %*% reg[nams,, drop = FALSE ]
  gg <- get_Z_X_XZ_formula(object$formulaTraits)
  traits0 <- stats::model.matrix(gg$formula_X0, constrasts = FALSE, data = object$data$dataTraits)
  msd = mean_sd_w(traits0, w = object$weights$columns)
  #mean_sd = cbind(mean = c(mean_sd1$mean),sd = c(mean_sd1$sd))
  #rownames(mean_sd) <- colnames(mean_sd1$mean)
  pred <- backscale_data(pred_scaled, msd)
  return(pred)

}

scale_data <- function(newdata, mean_sd){
  # newdata and mean_sd should  be matrices
  ones <- rep(1,nrow(newdata))
  Xc <- newdata - ones %*% t(mean_sd[,1,drop=FALSE][colnames(newdata),,drop= FALSE])
  Xc <- Xc / (ones %*% t(mean_sd[,2, drop =FALSE][colnames(newdata),,drop = FALSE]))
  return(Xc)
}

backscale_data <- function(pred_scaled, msd){
  # pred_scaled  should  be matrix and msd should be a data frame with a column mean and rownames
  ones <- rep(1,nrow(pred_scaled))
 # CWM <- CWM -  rep(1,nrow(CWM)) %*% msd$mean
  pred_scaled <- pred_scaled * (ones %*% msd$sd[1,colnames(pred_scaled),drop = FALSE])
  pred_scaled <- pred_scaled + ones %*% msd$mean[1,colnames(pred_scaled),drop = FALSE]
 return(pred_scaled)
}

predict_env <- function(object,  newdata, rank){
  # missing factors in newdata : no contribution of the missing factors
  reg <- predict_regr_traits(object, rank)
  nms <- rownames(reg)
  nms <- nms[nms %in% colnames(newdata)] # still needed for missing factors; see predict.dcca.r
  dat0 <- newdata[,nms, drop = FALSE]
  dat <- scale_data(dat0, mean_sd = object$c_traits_normed[,c(1,2)])
  #pred_scaled <- dat %*% reg[colnames(dat),, drop = FALSE ]
  # reg may not have the covariates (e.g. when dc-CA was started from CWM) and
  # newdata may not have all predictors
  nams <- intersect(colnames(dat), rownames(reg))
  pred_scaled <- dat[, nams,drop = FALSE] %*% reg[nams,, drop = FALSE ]
  gg <- get_Z_X_XZ_formula(object$formulaEnv)
  env0 <- stats::model.matrix(gg$formula_X0, constrasts = FALSE, data = object$data$dataEnv)
  msd = mean_sd_w(env0, w = object$weights$rows)
  #mean_sd = cbind(mean = c(mean_sd1$mean),sd = c(mean_sd1$sd))
  #rownames(mean_sd) <- colnames(mean_sd1$mean)
  pred <- backscale_data(pred_scaled, msd)
  return(pred)

}
