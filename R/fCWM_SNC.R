#' @title Calculate community weighted means and species niche centroids for 
#' double constrained correspondence analysis
#'
#' @description
#' Double constrained correspondence analysis (dc-CA) can be calculated directly 
#' from community weighted means (CWMs), with the trait data from which the 
#' CWMs are calculated, and the environmental data and weights for species 
#' and sites (the abundance totals for species and sites). Statistical testing
#' at the species level requires also the species niche centroids (SNCs).
#' The function \code{fCWM_SNC} calculates the CWMs and SNCs from the trait 
#' and environmental data, respectively, using a formula interface, so as to
#' allow categorical traits and environmental variables. The resulting object
#' can be set as the \code{response} argument in \code{\link{dc_CA}} so as to
#' give the same output as a call to \code{\link{dc_CA}} with the abundance
#' data as \code{response}, at least up to sign changes of the axes.
#'
#' @inheritParams dc_CA
#' 
#' @param minimal_output logical. Default \code{TRUE} for use of the return 
#' value as \code{response} in a call to \code{\link{dc_CA}}.
#' @param verbose logical for printing a simple summary (default: TRUE)
#
#' @details
#' The argument \code{formulaTraits} determines which CWMs are calculated.
#' The CWMs are calculated from the trait data (non-centered, non-standardized).
#' With trait covariates, the other predictor traits are adjusted for the trait 
#' covariates by weighted regression, after which the overall weighted mean 
#' trait is added. This has the advantage that each CWM has the scale of the 
#' original trait.
#'
#' The SNCs are calculated analogously from environmental data.
#'
#' Empty (all zero) rows and columns in \code{response} are removed from 
#' the \code{response} and the corresponding rows from \code{dataEnv} and 
#' \code{dataTraits}. Subsequently, any columns with missing values are 
#' removed from  \code{dataEnv} and \code{dataTraits}. It gives an error 
#' (object 'name_of_variable' not found), if variables with missing entries
#' are specified in \code{formulaEnv} and \code{formulaTraits}.
#'
#' In the current implementation, \code{formulaEnv} and \code{formulaTraits} 
#' should contain variable names as is, \emph{i.e.} transformations of 
#' variables in the formulas gives an error ('undefined columns selected') 
#' when the \code{\link{scores}} function is applied.
#'
#' @return
#' The default returns a list of CWM, SNC, weights, \code{formulaTraits} and
#' inertia (weighted variance explained by
#' the traits and by the environmental variables)
#' a list of data with elements \code{dataEnv} and \code{dataTraits}. When 
#' \code{minimal_output = FALSE}, some more statistics are given that are 
#' mainly technical or recomputed when the return value is used as 
#' \code{response} in a call to \code{\link{dc_CA}}.
#'
#' @references
#' Kleyer, M., Dray, S., Bello, F., Lepš, J., Pakeman, R.J., Strauss, 
#' B., Thuiller, W. & Lavorel, S. (2012) Assessing species and community 
#' functional responses to environmental gradients: which multivariate methods?
#' Journal of Vegetation Science, 23, 805-821.
#' \doi{10.1111/j.1654-1103.2012.01402.x}
#'
#' ter Braak, CJF, Šmilauer P, and Dray S. 2018. Algorithms and biplots for
#' double constrained correspondence analysis.
#' Environmental and Ecological Statistics, 25(2), 171-197.
#' \doi{10.1007/s10651-017-0395-x} 
#'
#' ter Braak C.J.F. and  P. Šmilauer  (2018). Canoco reference manual
#' and user's guide: software for ordination (version 5.1x).
#' Microcomputer Power, Ithaca, USA, 536 pp.
#'
#' Oksanen, J., et al. (2022)
#' vegan: Community Ecology Package. R package version 2.6-4.
#' \url{https://CRAN.R-project.org/package=vegan}.
#'
#' @seealso \code{\link{dc_CA}}, \code{\link{plot.dcca}}, 
#' \code{\link{scores.dcca}}, \code{\link{print.dcca}} and 
#' \code{\link{anova.dcca}}
#' 
#' @example demo/dune_fCWMSNC.r
#' @export
fCWM_SNC <- function(response = NULL, 
                     dataEnv = NULL, 
                     dataTraits = NULL,
                     formulaEnv = NULL,
                     formulaTraits = NULL,
                     divideBySiteTotals = TRUE,
                     dc_CA_object = NULL,
                     minimal_output = TRUE,
                     verbose = TRUE) {
  # response matrix or data frame, dataEnv and dataTraits data frames in which 
  # formualaE and formulaT are evaluated
  # dc_CA_object = result (value) of a previous run, can be used to save 
  # computing time for runs that modify the formula for samples 
  # (step2: RDAonEnv) only
  # The step1 (CCAonTraits and the data and formulaTraits) are taken from 
  # dc_CA_object into the new result.
  # If set, formulaTraits, response, dataEnv, dataTraits are not used at all 
  # and have no efffect on the result
  call <- match.call()
  out <- check_data_dc_CA(formulaEnv, formulaTraits,
                          response, dataEnv, dataTraits, divideBySiteTotals, call)
  # CWM and CWM ortho
  ms <- try(f_wmean(out$formulaTraits, tY = out$data$Y, out$data$dataTraits,
                weights=out$weights, name= "CWM"))
  if (inherits(ms, "try-error")) {
    stop("singular trait data. No CWMs generated.\n") 
  }   
  # SNC and SNC ortho
  mt <- try(f_wmean(out$formulaEnv, tY = t(out$data$Y), out$data$dataEnv,
                weights=out$weights, name= "SNC"))
  if (inherits(mt, "try-error")) {
    stop("singular environment data. No SNCs generated\n") 
  } 
  out <- list(
    CWM = ms$wmean,
    SNC = mt$wmean,
    formulaEnv = out$formulaEnv,
    formulaTraits = out$formulaTraits, 
    inertia = c(traits_explain=ms$explained, env_explain =mt$explained),
    weights = out$weights, 
    call = out$call, 
    data = out$data[-1]# remove Y from out$data
    )
if (!minimal_output) { 
     out <- c(out, list( 
      CWMs_orthonormal_traits = ms$wmean_ortho,
      SNCs_orthonormal_env = mt$wmean_ortho,
      trans2ortho = list(CWM2CWM_ortho = ms$to_ortho, 
                         SNC2SNC_ortho = mt$to_ortho),
      T_ortho = ms$to_ortho,
      E_ortho = mt$to_ortho
    ))
}
  return(out)
}

#' @noRd
#' @keywords internal
f_canonical_coef_traits2 <- function(out){
  # do first the trivia of mean, sds and VIF
  # analogously  to calculate_b_se_tval in f_dc_CA.r
  XZ <- is.null(out$CWM2CWM_ortho)
  ms <- msdvif(formula = out$formulaTraits, data = out$data$dataTraits, 
               weights = out$weights$columns, XZ = XZ)
  avg <- ms$meansdvif[, 1]
  sds <- ms$meansdvif[, 2]
  VIF <- ms$meansdvif[, 3]
  # end of trivia of mean sds and VIF
  if (!is.null(out$CWM2CWM_ortho)) {
    B1_traitsN <- out$CWM2CWM_ortho * sds
  } else {
    if (inherits(out$CCAonTraits, "cca")) {
      B1_traitsN <- scores(out$CCAonTraits, display = "reg", 
                           scaling = "species", 
                           choices = 1:rank_mod(out$CCAonTraits))
    } else {
      return(ms$meansdvif)
    }
  }
  if (inherits(out$RDAonEnv, "wrda")) {
    B_star <- scores(out$RDAonEnv, display = "species", scaling = "sites", 
                     choices = 1:rank_mod(out$RDAonEnv))
  } else  {
    B_star <- scores(out$RDAonEnv, display = "species", 
                     scaling = "sites", choices = 1:rank_mod(out$RDAonEnv), 
                     const = 1)
  }
  # End appendix 6.2 in ter Braak  et al 2018
  c_traitsN <- B1_traitsN %*% B_star
  colnames(c_traitsN) <- paste0("Regr", seq_len(ncol(c_traitsN)))
  coef_normed <- cbind(Avg = avg, SDS = sds, VIF = VIF, c_traitsN, tval1 = NA)
  if (is.na(avg[1])) {
    coef_normed[, "SDS"] <- NA
    attr(coef_normed, "warning") <-
      paste("regression coefs can only be interpreted columnwise in terms of ", 
            "their sign and quantitatively only rowwise ", 
            "as the standard deviation of traits is NA." )
  } else {
    attr(coef_normed, "meaning") <- 
      "mean, sd, VIF, standardized regression coefficients, t-values not available"
  }
  return(coef_normed)
}

#' @noRd
#' @keywords internal
f2_orth <- function(CWM,
                    formulaTraits, 
                    dataTraits, 
                    weights.cols, 
                    weights.rows,
                    name = "CWM") {
  msqr <- msdvif(formulaTraits, dataTraits, weights.cols, XZ = FALSE)
  sWn <- sqrt(weights.cols)
  X <- msqr$Xw/sWn 
  msd <- msqr$meansdvif
  CWM2CWM_ortho <- solve(qr.R(msqr$qrX))
  colnames(CWM2CWM_ortho) <- rownames(CWM2CWM_ortho)
  if (all(rownames(CWM2CWM_ortho) %in% colnames(CWM))) {
    CWM <- as.matrix(CWM[, rownames(CWM2CWM_ortho), drop = FALSE])
  } else {
    if (name == "CWM") {
      fname <- "formulaTraits"
    } else {
      fname <- "formulaEnv"
    }
    stop("All names generated by ", fname, " namely \n",
         paste(rownames(CWM2CWM_ortho), collapse = ","), 
         "\n should occur in response$", name, " being:\n", 
         paste(colnames(CWM), collapse = ","), "\n" )
  }
  msd <- mean_sd_w(CWM, w = weights.rows)
  CWM <- CWM - rep(1, nrow(CWM)) %*% msd$mean
  CWMs_orthonormal_traits <- 
    CWM[, rownames(CWM2CWM_ortho), drop = FALSE] %*% CWM2CWM_ortho
  rownames(CWMs_orthonormal_traits) <- rownames(CWM)
  return(list(CWMs_orthonormal_traits = CWMs_orthonormal_traits, 
              CWM2CWM_ortho = CWM2CWM_ortho))
}

#' @noRd
#' @keywords internal
checkCWM2dc_CA <- function(object, 
                           dataEnv,
                           dataTraits,
                           formulaTraits) {
  # object is from CWMSNC object
  if (!is.null(formulaTraits)) {
    object$formulaTraits <- formulaTraits
  } else if (is.null(object$formulaTraits)) {
    warning("formulaTraits set to ~. in checkCWM2dc_CA.\n")
    object$formulaTraits <- ~.
  }
  object$CWM <- as.data.frame(object$CWM)
  if (!is.null(object$response$SNC)) {
    object$response$SNC <- as.data.frame(object$response$SNC)
  }
  object$Nobs <- nrow(object$CWM)
  # check object
  ## check data dataEnv dataTraits
  if (is.null(object[["data"]])) {
    object[["data"]] <- list()
  }
  if ("dataTraits" %in% names(object)) {
    if (is.null(dataTraits)) {
      dataTraits<- as.data.frame(object$dataTraits)
      dataTraits <- as.data.frame(lapply(dataTraits, function(x) {
        if (is.character(x)) x<- as.factor(x) else x
        return(x)
      }))
      rownames(dataTraits) <- rownames(object$dataTraits)
      object[["data"]]$dataTraits <- dataTraits
      
    }
    object$dataTraits <- NULL
  }
  if ("dataEnv" %in% names(object)) {
    if (is.null(dataEnv)) {
      dataEnv<- as.data.frame(object$dataEnv)
      dataEnv <- as.data.frame(lapply(dataEnv, function(x) {
        if (is.character(x)) x<- as.factor(x) else x
        return(x) 
      }))
      rownames(dataEnv) <- rownames(object$dataEnv)
      object[["data"]]$dataEnv <- object$dataEnv
    }
    object$dataEnv<-NULL
  }
  if (is.null(dataEnv)) {
    if (is.null(object$data$dataEnv)) {
      stop("Supply environmental data to the dc_CA function.\n")
    }
  } else {
    object$data$dataEnv <- as.data.frame(lapply(dataEnv, function(x){
      if (is.character(x)) x <- as.factor(x) else x
      return(x) 
    }))
  }
  if (is.null(dataTraits)) {
    if (!is.null(object$dataTraits)) {
      object$data$dataTraits <- object$dataTraits 
    } else if (is.null(object$data$dataTraits)) {
      warning("Supply trait data to the dc_CA function.\n")
    }
  } else { 
    warning("Trait data taken from the argument dataTraits.\n")
    object$data$dataTraits <- as.data.frame(lapply(dataTraits, function(x) {
      if (is.character(x)) x <- as.factor(x) else x
      return(x) 
    }))
  }
  # check weights
  if (is.null(object$weights)) {
    # try dataTraits and dataEnv
    if (!is.null(object$data$dataTraits$weight)) {
      warning("species weights taken from dataTraits$weight.\n")
      object$weights <- c(list(columns = object$data$dataTraits$weight),
                          object$weights)
    }
    if (!is.null(object$data$dataEnv$weight)) {
      warning("site weights taken from dataEnv$weight.\n")
      object$weights <- c(object$weights, 
                          list(rows = object$data$dataEnv$weight))
    }
  }
  if (is.null(object$weights)) {
    warning("no weights supplied with response$CWM; weigths all set to 1.\n")
    object$weights <- list(columns = rep(1 / nrow(object$data$dataTraits),
                                         nrow(object$data$dataTraits)),
                           rows = rep(1 / object$Nobs, object$Nobs)
                           )
  } else if (!is.list(object$weights)) {
    ll <- length(object$weights)
    if (ll == object$Nobs) {
      warning("no species weights supplied with response$CWM; ",
              "weigths all set to 1.\n")
      object$weights <- list(columns = object$weights,
                        rows = rep(1 / object$Nobs, object$Nobs)
                             )
    } else if (ll == nrow(object$data$dataTraits)) {
      warning("no site weights supplied with response$CWM; ", 
              "site weigths all set to 1.\n")
      object$weights <- list(
        columns = rep(1 / nrow(object$data$dataTraits),
                                          nrow(object$data$dataTraits)),
        rows = object$weights
                             )
    }
  }
  if (!is.null(object$weights)){
    if (is.null(object$weights$rows)){
      warning("no site weights supplied with response$CWM; ", 
              "site weigths all set to 1.\n")
      object$weights$rows <- rep(1 / object$Nobs, object$Nobs)
    }
    if (is.null(object$weights$columns)) {
      warning("no species weights supplied with response$CWM; ",
              "species weigths all set to 1.\n")
      object$weights$columns <- rep(1 / nrow(object$data$dataTraits),
                                    nrow(object$data$dataTraits))
    }
  }
  # make sure columns is the first in the list...
  object$weights <- list(columns = object$weights$columns, rows = object$weights$rows)
  object$weights <- lapply(object$weights, function(x)x/sum(x))
  # change ~. to names
  formulaTraits <- change_reponse(object$formulaTraits, "Y", 
                                  object$data$dataTraits)
  # object checked..
  CWM2ortho <- with(object, f2_orth(CWM, formulaTraits, data$dataTraits, 
                                    weights$columns, weights$rows, 
                                    name = "CWM"))
  object$CWMs_orthonormal_traits <- 
    CWM2ortho$CWMs_orthonormal_traits * sqrt((object$Nobs - 1) / object$Nobs)
  if (!is.null(object$SNC) && !is.null(object$weights$rows)) {
    SNC2ortho <- with(object, f2_orth(SNC, formulaEnv, data$dataEnv, 
                                      weights$rows, weights$columns, 
                                      name = "SNC"))
    object$SNCs_orthonormal_env <- SNC2ortho$CWMs_orthonormal_traits
  }
  object$CWM2CWM_ortho <- CWM2ortho$CWM2CWM_ortho
  return(object)
}
#' @noRd
#' @keywords internal
f_wmean <- function(formulaEnv, tY, dataEnv, weights, name = "SNC"){
  # SNC and SNC ortho or 
  # CWM and CWM ortho if name = "CwM" from formulaTraits, Y and dataTraits
  tot <- sum(tY)
  formulaEnv <- change_reponse(formulaEnv, "Y", dataEnv)
  if (name == "SNC")  w <- weights$rows else w <- weights$columns
  msqr <- msdvif(formulaEnv, dataEnv, w, XZ = FALSE)
  sWn <- sqrt(w)
  X <- msqr$Xw / sWn 
  msd <- msqr$meansdvif
  E_ortho <- qr.Q(msqr$qrX) / sWn
  # so Q represents the orthonormalized predictors
  if (name == "SNC")  w <- weights$columns else w <- weights$rows
  SNC <- diag(1 /( tot * w)) %*% as.matrix(tY) %*% X
  SNC2SNC_ortho <- try(solve(qr.R(msqr$qrX)))
  if (inherits(SNC2SNC_ortho, "try-error")) {
    warning("singular environment data. Env_explain not available.\n") 
    SNC2SNC_ortho <- matrix(NA,nrow = ncol(SNC), ncol= 1)
  }
  SNCs_orthonormal_env <- SNC %*% SNC2SNC_ortho
  env_explain <- sum(SNCs_orthonormal_env ^ 2 * w)
  SNC <- SNC + rep(1, nrow(SNC)) %*% matrix(msd[, 1], nrow = 1)
  rownames(SNC) <- rownames(SNCs_orthonormal_env) <- rownames(tY)
  rownames(E_ortho) <- colnames(tY)
  res <- list(wmean = SNC, explained = env_explain,
              wmean_ortho = SNCs_orthonormal_env, name = name, to_ortho=E_ortho)
  return(res)
}
#' @noRd
#' @keywords internal
check_data_dc_CA <- function(formulaEnv, formulaTraits,
                             response, dataEnv, dataTraits, divideBySiteTotals, call){
  #  check and amend: make sure there are no empty rows or columns
  if (any(is.na(response))) {
    stop("The response should not have missing entries.\n")
  }
  if (any(response < 0)) {
    stop("The response should not have negative values.\n")
  }
  if (is.null(dataEnv)) {
    stop("dataEnv must be specified in dc_CA.\n")
  } else  dataEnv <- as.data.frame(dataEnv)
  if (is.null(dataTraits)) {
    stop("dataTraits must be specified in dc_CA.\n")
  } else  dataTraits <- as.data.frame(dataTraits)
  if (!is.matrix(response)) {
    response <- as.matrix(response)
  }
  id0 <- 1
  while(length(id0)) {
    TotR <- rowSums(response)
    id0 <- which(TotR == 0)
    if (length(id0)) {
      response <- response[-id0, ]
      dataEnv  <- dataEnv[-id0, ]
    }
    TotC <- colSums(response)
    id0 <- which(TotC == 0)
    if (length(id0)) {
      response <- response[, -id0]
      dataTraits <- dataTraits[-id0, ]
    }
  }
  # delete columns with missing data
  id <- rep(FALSE, ncol(dataEnv))
  for (ii  in seq_along(id)) {
    id[ii] <- sum(is.na(dataEnv[, ii])) == 0
    if (!id[[ii]]) {
      warning("variable", names(dataEnv)[ii],
              "has missing values and is deleted from the environmental data.\n")
    }
  }
  dataEnv <- dataEnv[, id]
  id <- rep(FALSE, ncol(dataTraits))
  for (ii  in seq_along(id)) {
    id[ii] <- sum(is.na(dataTraits[,ii])) == 0
    if (!id[[ii]]) {
      warning("variable", names(dataTraits)[ii],
              "has missing values and is deleted from trait data.\n")
    }
  }
  dataTraits <- dataTraits[, id]
  dataEnv <- as.data.frame(lapply(dataEnv, function(x) {
    if (is.character(x)) x <- as.factor(x) else x
    return(x)
  }))
  dataTraits <- as.data.frame(lapply(dataTraits, function(x) {
    if (is.character(x)) x <- as.factor(x) else x
    return(x)
  }))
  rownames(dataEnv) <- rownames(response)
  rownames(dataTraits) <- colnames(response)
  # end of check
  TotR <- rowSums(response)
  if (divideBySiteTotals) {
    response <- response / (TotR %*% t(rep(1, ncol(response))))
    TotR <- rep(1, length(TotR))
  }
  TotC <- colSums(response)
  TotR <- TotR /sum(TotR)
  TotC <- TotC/sum(TotC)
  weights <- list(columns = TotC, rows = TotR) # unit sums
  Nobs <- nrow(response)
  if (is.null(formulaTraits)) {
    formulaTraits <- as.formula(paste("~", paste0(names(dataTraits),
                                                  collapse = "+")))
    warning("formulaTraits set to ~. in fCWMSNC.\n")
  }
  if (is.null(formulaEnv)) {
    formulaEnv <- as.formula(paste("~", paste0(names(dataEnv),
                                               collapse = "+")))
    warning("formulaEnv set to ~. in fCWMSNC.\n")
  }
  
  out <- list(formulaTraits=formulaTraits,formulaEnv= formulaEnv,
              data = list(Y=response, dataEnv= dataEnv, dataTraits=dataTraits),
              call = call,  weights = weights, Nobs = Nobs)
  return(out)
}
