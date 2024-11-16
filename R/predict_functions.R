#' @noRd
#' @keywords internal
predict_regr_env <- function(object, 
                             rank, normed = TRUE) {
  sc <- scores_dcca(object, choices = seq_len(rank),
          display = c("reg", "bp_traits"), normed = normed)
  B_env_regr <- sc[["regression"]][, -c(1, 2, 3), drop = FALSE]
  C_trait_bip <- sc$biplot_traits
  if (normed)standardized <- "standardized " else standardized <- ""
  reg <- B_env_regr %*% t(C_trait_bip)
  attr(reg, which = "meaning") <-paste0(standardized,
    "regression coefficients to predict traits from environment.")
  if (!normed) 
    attr(reg, which = "mean")<- list(mean_traits = attr(sc$biplot_traits, which = "mean"),
                                     mean_env = sc[["regression"]][,"Avg"])  
  return(reg)
}

#' @noRd
#' @keywords internal
predict_regr_traits <- function(object, 
                                rank, normed = TRUE) {
  sc <- scores_dcca(object, choices = seq_len(rank), 
            display = c("reg_traits", "bp"), normed = normed)
  B_traits_regr <- sc$regression_traits[, -c(1, 2, 3), drop = FALSE]
  C_env_bip <- sc[["biplot"]]
  if (normed) standardized <- "standardized " else standardized <- ""
  reg <- B_traits_regr %*% t(C_env_bip) 
  attr(reg, which = "meaning") <-paste0(standardized,
          "regression coefficients to predict environment from traits.")
  if (!normed) 
    attr(reg, which = "mean")<- list(mean_traits = sc$regression_traits[,"Avg"],
                                     mean_env = attr(sc[["biplot"]], which = "mean")) 
  return(reg)
}

#' @noRd
#' @keywords internal
predict_traits <- function(object, 
                           newdata1, 
                           rank) {
  # missing factors in newdata1 : the reference level (first level) is assumed
  reg <- predict_regr_env(object, rank, normed = FALSE)
  reg[is.na(reg)] <- 0
  #pred_mean <- fpred_scaled(newdata1, reg)
  newdata1 <- set_newdata( object,newdata1,type = "traitsFromEnv",  
                           means_mis = attr(reg, which = "mean")$mean_env)
  # NB: newdata1 derived from model.matrix!
  newdata1 <- subtract_mean(newdata1, mean0 =attr(reg, which = "mean")$mean_env)
  pred_mean <- newdata1 %*% reg
  pred <- add_mean(pred_mean, mean0 =attr(reg, which = "mean")$mean_traits)
  return(pred)
}

#' @noRd
#' @keywords internal
subtract_mean <- function(dat0, mean0) {
 # dat0 and mean0 should be a vector means(substract mean only)
  ones <- rep(1, nrow(dat0))
  Xc <- dat0 - ones %*% t(mean0)
  return(Xc)
}

#' @noRd
#' @keywords internal
add_mean <- function(pred_mean,
                           mean0) {
  # pred_mean  should  be matrix and mean0 should be a vector means
  
  ones <- rep(1, nrow(pred_mean))
 # pred_mean <- pred_mean * 
 #   (ones %*% msd$sd[1, colnames(pred_mean), drop = FALSE])
  pred_mean <- pred_mean + ones %*% t(mean0)
  return(pred_mean)
}

#' @noRd
#' @keywords internal
set_newdata <- function(object,newdata, type, means_mis) {
  # check for 1 data frame (either env or traits)
  # where check_newdata gives environmental data and trait data respectively
  if (type %in%  c("traitsFromEnv", "CWM")) {
  # print(names(object$data$dataEnv))
  #  print(object$formulaEnv)

    ff_get <- get_Z_X_XZ_formula(object$formulaEnv, object$data$dataEnv)
    #if (is.null(newdata)) newdata1 <- object$data$dataEnv else newdata1 <- newdata
    trainingData <- object$data$dataEnv
  } else if (type %in% c("envFromTraits", "SNC")) { # "env", "reg_traits"
    ff_get <- get_Z_X_XZ_formula(object$formulaTraits, object$data$dataTraits)
    #if (is.null(newdata))   newdata1 <- object$data$dataTraits else newdata1 <-newdata
    trainingData <- object$data$dataTraits
  }
  missing_names <- !ff_get$all_nams %in% names(newdata)
  if (any(missing_names)){
    missing_names <- ff_get$all_nams[missing_names]
    warning("newdata does not contain the predictor variables\n ", 
            paste(missing_names, collapse = ","),
            "\nThese are set at their mean values and,\n", 
            "for factors, at the reference level",
            paste("\nThe current formula is\n",
            as.character(ff_get$formula_XZ )[1],
            as.character(ff_get$formula_XZ )[2]) )
            
    # add variables /factors to newdata1
    newdf <- data.frame(matrix(NA,nrow = nrow(newdata), ncol = length(missing_names),
                               dimnames=list(rownames(newdata), missing_names)))
    
    for (nam in missing_names){
      # for quantitative factors set value to the mean
       if (nam %in% names(means_mis)) newdf[[nam]] <-means_mis[nam]
      }
    newdata <- cbind(newdata, newdf)
    
  } # end if any mising
  
  for (nam in c(ff_get$focal_factor, ff_get$Condi_factor) ){
    newdata[[nam]] <- factor(newdata[[nam]],levels = levels(trainingData[[nam]]));
    newdata[[nam]][is.na(newdata[[nam]])] <- levels(trainingData[[nam]])[1]
  }

  newdata <- model.matrix(ff_get$formula_XZ, 
                       data = newdata)[, -1, drop = FALSE]
  if (!nrow(newdata)) warning("newdata does not contain informative data.")
  return(newdata)
}

#' @noRd
#' @keywords internal
predict_env <- function(object,
                        newdata1, 
                        rank) {
  # missing factors in newdata1 : the reference level (first level) is assumed
  reg <- predict_regr_traits(object, rank, normed = FALSE)
  reg[is.na(reg)] <- 0
  #pred_mean <- fpred_scaled(newdata1, reg)
  newdata2 <- set_newdata(object,newdata1,  type = "envFromTraits", 
                          means_mis = attr(reg, which = "mean")$mean_traits)
  newdata1 <- subtract_mean(newdata2, mean0 =attr(reg, which = "mean")$mean_traits)
  pred_mean <- newdata1 %*% reg
  pred <- add_mean(pred_mean, mean0 =attr(reg, which = "mean")$mean_env)
  return(pred)
}

#' @noRd
#' @keywords internal
predict_response <- function(object, 
                             newdata1,
                             rank, weights = object$weights) {
  # newdata1 must be a list two dataframes, element 1: trait and  element 2 env data
  sc <- scores_dcca(object, 
        choices = seq_len(rank), display = c("reg", "reg_traits"), scaling = "sym",
        normed = FALSE)
  B_traits_regr <- sc[["regression_traits"]][, -c(1, 2, 3), drop = FALSE] 
  B_traits_regr[is.na(B_traits_regr)]<-0
  newdata2 <- set_newdata( object,newdata1[[1]],type = "envFromTraits", 
                           means_mis = sc[["regression_traits"]][,"Avg"])
  pred_scaled_species <- subtract_mean(newdata2,
                                        mean0 = sc[["regression_traits"]][,"Avg"])
  pred_scaled_species <- 
    pred_scaled_species %*% B_traits_regr
  B_env_regr <- sc[["regression"]][, -c(1, 2, 3), drop = FALSE]
  B_env_regr[is.na(B_env_regr)]<-0
  # attr(reg, which = "mean")<- list(mean_traits = attr(sc$biplot_traits, which = "mean"),
  #                                  mean_env = sc[["regression"]][,"Avg"])  
  newdata2 <- set_newdata( object,newdata1[[2]],type = "traitsFromEnv", 
                           means_mis = sc[["regression"]][,"Avg"])
  pred_scaled_sites <- subtract_mean(newdata2, 
                                        mean0 = sc[["regression"]][,"Avg"])
  pred_scaled_sites <- pred_scaled_sites%*% B_env_regr
  interact <- pred_scaled_sites %*% t(pred_scaled_species)
  pred <- (1 + interact) * 
    ( rep(1,nrow(interact)) %*% t(weights[[1]]/sum(weights[[1]])) )* weights[[2]]
  return(pred)
}
#' @noRd
#' @keywords internal
predict_regr_all <- function(object, 
                             rank) {
  # regresion coefficient of transformed response on env and trait predictors
  # value is env by trait
  sc <- scores(object, choices = seq_len(rank), display = c("reg", "reg_traits"))
  reg <-sc[["regression_traits"]][,-c(1,2,3), drop = FALSE]%*% 
    t(sc[["regression"]][,-c(1,2,3), drop = FALSE])
  attr(reg, which = "meaning") <-
  "regression coefficients for traits and environment to predict the response."
  return(reg)
}
#' @noRd
#' @keywords internal
predict_fc <- function(object, 
                                   rank) {
  # fourth-corner coefficients of transformed response on env and trait predictors
  # value is env by trait
  sc <- scores(object, choices = seq_len(rank), display = c("bp", "bp_traits"))
  fc <- sc[["biplot_traits"]] %*% t(sc[["biplot"]])
  attr(fc, which = "meaning") <- "fourth-corner correlation"
  return(fc)
}
