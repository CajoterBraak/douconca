#' @noRd
#' @keywords internal
predict_regr_env <- function(object, 
                             rank) {
  sc <- scores(object, choices = seq_len(rank), display = c("reg", "bp_traits"))
  B_env_regr <- sc$regression
  C_trait_bip <- sc$biplot_traits
  reg <- B_env_regr[, -c(1, 2, 3), drop = FALSE] %*% t(C_trait_bip)
  return(reg)
}

#' @noRd
#' @keywords internal
predict_regr_traits <- function(object, 
                                rank) {
  sc <- scores(object, choices = seq_len(rank), display = c("reg_traits", "bp"))
  B_traits_regr <-sc$regression_traits
  C_env_bip <- sc$biplot
  reg <- B_traits_regr[, -c(1, 2, 3), drop = FALSE] %*% t(C_env_bip)
  return(reg)
}

#' @noRd
#' @keywords internal
predict_traits <- function(object, 
                           newdata1, 
                           rank) {
  # missing factors in newdata1 : no contribution of the missing factors
  reg <- predict_regr_env(object, rank)
  reg[is.na(reg)] <- 0
  pred_scaled <- fpred_scaled(newdata1, reg)
  gg <- get_Z_X_XZ_formula(object$formulaTraits)
  traits0 <- modelmatrixI(formula = gg$formula_X1, 
                          data = object$data$dataTraits, XZ = FALSE)
  msd <- mean_sd_w(traits0, w = object$weights$columns)
  pred <- backscale_data(pred_scaled, msd)
  return(pred)
}

#' @noRd
#' @keywords internal
scale_data <- function(dat0,
                       mean_sd) {
  # dat0 and mean_sd should  be matrices
  nams <- intersect(colnames(dat0), rownames(mean_sd))
  ones <- rep(1, nrow(dat0))
  Xc <- dat0[, nams, drop = FALSE] - ones %*% t(mean_sd[nams, 1, drop = FALSE])
  Xc <- Xc / (ones %*% t(mean_sd[nams, 2, drop = FALSE]))
  return(Xc)
}

#' @noRd
#' @keywords internal
backscale_data <- function(pred_scaled,
                           msd) {
  # pred_scaled  should  be matrix and msd should be a data frame with a column mean and rownames
  ones <- rep(1, nrow(pred_scaled))
  pred_scaled <- pred_scaled * 
    (ones %*% msd$sd[1, colnames(pred_scaled), drop = FALSE])
  pred_scaled <- pred_scaled + 
    ones %*% msd$mean[1, colnames(pred_scaled), drop = FALSE]
  return(pred_scaled)
}

#' @noRd
#' @keywords internal
check_newdata <- function(object, 
                          newdata,
                          type) {
  # check for 1 data frame (either env or traits)
  # BEWARE,
  # type is "traits" or "env"
  # where check_newdata gives environmental data and trait data respectively
  if (type == "traits") {
    wmff <- msdvif(object$formulaEnv, object$data$dataEnv, object$weights$rows,
                   XZ = TRUE, novif = TRUE)
    wm <- wmff$msd$mean
    ff <- wmff$ff_get
    c_normed <- cbind(Avg = c(wmff$msd$mean), SDS = c(wmff$msd$sd))
    rownames(c_normed) <- colnames(wmff$msd$mean)
    if (is.null(newdata)) newdata1 <- object$data$dataEnv
  } else if (type == "env") { # "env", "reg_traits"
    ff <- get_Z_X_XZ_formula(object$formulaTraits)
    wm <- t(object$c_traits_normed0[, "Avg", drop = FALSE])
    c_normed <- object$c_traits_normed0
    if (is.null(newdata)) {
      newdata1 <- object$data$dataTraits
    }
  }
  if (!is.null(newdata)) {
    newdata1 <- newdata
    nams <- !ff$all_nams %in% names(newdata1)
    nams <- ff$all_nams[nams]
    for (n in nams) {
      if (n %in% colnames(wm)) {
        newdata1[[n]] <- wm[1,n] # set to the mean
      } else {
        newdata1[[n]] <- 0 
      }
    }
  }
  dat0 <- model.matrix(ff$formula_XZ, constrasts = FALSE, 
                       data = newdata1)[, -1, drop = FALSE]
  newdata1 <- scale_data(dat0, mean_sd = c_normed[, c(1, 2), drop = FALSE])
  return(newdata1)
}

#' @noRd
#' @keywords internal
predict_env <- function(object,
                        newdata1, 
                        rank) {
  # missing factors in newdata1 : no contribution of the missing factors
  reg <- predict_regr_traits(object, rank)
  reg[is.na(reg)] <- 0
  pred_scaled <- fpred_scaled(newdata1, reg)
  gg <- get_Z_X_XZ_formula(object$formulaEnv)
  env0 <- modelmatrixI(formula = gg$formula_X1, data = object$data$dataEnv, 
                       XZ = FALSE)
  msd <- mean_sd_w(env0, w = object$weights$rows)
  pred <- backscale_data(pred_scaled, msd)
  return(pred)
}

#' @noRd
#' @keywords internal
fpred_scaled <- function(newdata1, 
                         reg) {
  nms <- rownames(reg)
  # still needed for missing factors; see predict.dcca.r
  nms <- nms[nms %in% colnames(newdata1)] 
  # reg may not have the covariates (e.g. when dc-CA was started from CWM) and
  # newdata1 may not have all predictors
  nams <- intersect(colnames(newdata1), rownames(reg))
  reg[is.na(reg)] <- 0
  pred_scaled <- newdata1[, nams, drop = FALSE] %*% reg[nams, , drop = FALSE]
  return(pred_scaled)
}

#' @noRd
#' @keywords internal
predict_response <- function(object, 
                             newdata1,
                             rank) {
  # newdata1 must be a list two dataframes, element 1: trait and  element 2 env data
  B_traits_regr <- scores(object, choices = seq_len(rank), 
                          display = "reg_traits", scaling  = "sites")
  pred_scaled_species <- 
    fpred_scaled(newdata1[[1]], B_traits_regr[, -c(1, 2, 3), drop = FALSE])
  B_env_regr <- scores(object, choices = seq_len(rank), display = "reg", 
                       scaling = "sites")
  pred_scaled_sites <- fpred_scaled(newdata1[[2]], 
                                    B_env_regr[, -c(1, 2, 3), drop = FALSE])
  interact <- pred_scaled_sites %*% t(pred_scaled_species)
  pred <- (1 + interact) * 
    rep(1, rep(nrow(interact))) %*% t(object$weights$columns)
  return(pred)
}
