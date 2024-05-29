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

predict_traits <- function(object, newdata1, rank){
  # missing factors in newdata1 : no contribution of the missing factors
  reg <- predict_regr_env(object, rank)
  pred_scaled <- fpred_scaled(newdata1, reg)
  gg <- get_Z_X_XZ_formula(object$formulaTraits)
  traits0 <- stats::model.matrix(gg$formula_X0, data = object$data$dataTraits)
  msd = mean_sd_w(traits0, w = object$weights$columns)
  #mean_sd = cbind(mean = c(mean_sd1$mean),sd = c(mean_sd1$sd))
  #rownames(mean_sd) <- colnames(mean_sd1$mean)
  pred <- backscale_data(pred_scaled, msd)
  return(pred)

}

scale_data <- function(newdata1, mean_sd){
  # newdata1 and mean_sd should  be matrices
  ones <- rep(1,nrow(newdata1))
  Xc <- newdata1 - ones %*% t(mean_sd[,1,drop=FALSE][colnames(newdata1),,drop= FALSE])
  Xc <- Xc / (ones %*% t(mean_sd[,2, drop =FALSE][colnames(newdata1),,drop = FALSE]))
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

check_newdata <- function(object, newdata, type){
  # check for 1 data frame (either env or traits)
  # BEWARE,
  # type is "traits" or "env"
  # where check_newdata gives environmental data and trait data respectively

  if (type %in% c("traits")) {
    #ff <- get_Z_X_XZ_formula(object$formulaEnv)
    wmff<- msdvif(object$formulaEnv, object$data$dataEnv,  object$weights$rows,
                XZ = TRUE, novif= TRUE)
    wm <- wmff$mean; ff <- wmff$ff_get
    c_normed <- cbind(Avg = c(wmff$msd$mean), SDS = c(wmff$msd$sd) )
    rownames(c_normed) <- colnames(wmff$msd$mean)
    if (is.null(newdata)) newdata1 <- object$data$dataEnv
  } else if (type %in% c("env")){ # "env", "reg_traits"
    ff <- get_Z_X_XZ_formula(object$formulaTraits)
    wm <- t(object$c_traits_normed0[,"Avg", drop = FALSE])
    c_normed <- object$c_traits_normed0
    if (is.null(newdata)) newdata1 <- object$data$dataTraits
  }
  if (!is.null(newdata)){
    newdata1 <- newdata
    nams <- !colnames(wm) %in% names(newdata1)
    nams <-colnames(wm)[nams]
    for (n in nams){
      if (n %in% colnames(wm)) newdata1[[n]] <- wm[1,n] # set to the mean
          #c_normed[n,"Avg"] # set to average
      else newdata1[[n]] <- 0 # for factors
    }

  }
  dat0 <- stats::model.matrix(ff$formula_XZ, constrasts = FALSE, data = newdata1)[,-1,drop=FALSE]
  #colnames(dat0)
  newdata1  <- scale_data(dat0, mean_sd = c_normed[,c(1,2)])
  return(newdata1)
}

predict_env <- function(object,  newdata1, rank){
  # missing factors in newdata1 : no contribution of the missing factors
  reg <- predict_regr_traits(object, rank)
  pred_scaled <- fpred_scaled(newdata1, reg)
  gg <- get_Z_X_XZ_formula(object$formulaEnv)
  env0 <- stats::model.matrix(gg$formula_X0, data = object$data$dataEnv)
  msd = mean_sd_w(env0, w = object$weights$rows)
  #mean_sd = cbind(mean = c(mean_sd1$mean),sd = c(mean_sd1$sd))
  #rownames(mean_sd) <- colnames(mean_sd1$mean)
  pred <- backscale_data(pred_scaled, msd)
  return(pred)

}

fpred_scaled <- function(newdata1, reg){
  #reg <- predict_regr_traits(object, rank)
  nms <- rownames(reg)
  nms <- nms[nms %in% colnames(newdata1)] # still needed for missing factors; see predict.dcca.r
  # reg may not have the covariates (e.g. when dc-CA was started from CWM) and
  # newdata1 may not have all predictors
  nams <- intersect(colnames(newdata1), rownames(reg))
  pred_scaled <- newdata1[, nams,drop = FALSE] %*% reg[nams,, drop = FALSE ]
  return(pred_scaled)
}

predict_response <- function(object, newdata1, rank){
  # newdata1 must be a list two dataframes, element 1: trait and  element 2 env data
  B_traits_regr <- scores(object, choices = seq_len(rank), display = c("reg_traits"), scaling ="sites")
  pred_scaled_species <- fpred_scaled(newdata1[[1]], B_traits_regr[,-c(1,2,3), drop = FALSE])
  B_env_regr <- scores(object, choices = seq_len(rank), display = c("reg"), scaling ="sites")
  pred_scaled_sites <- fpred_scaled(newdata1[[2]], B_env_regr[,-c(1,2,3), drop = FALSE])
  interact <- pred_scaled_sites %*% t(pred_scaled_species)
  pred <- (1+interact) * rep(1, rep(nrow(interact)))%*% t(object$weights$columns)
  return(pred)
}


