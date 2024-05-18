#' @title Calculate community weighted means and species niche centroids for double constrained correspondence analysis
#'
#' @description
#' Double constained correspondence analysis (dc-CA) can be calculated directly from
#' community weighted means (CWMs), with the trait data from
#' which the CWMs are calculated, and the enviromental data.
#' Statistical testing at the species level requires also the species niche centroids (SNCs).
#' The function \code{fCWM_SNC} calculates the CWMs and SNCs from the trait and environmental
#' data, respectively, using a formula interface, so as to allow categorical traits and
#' environmental variables.
#' The resulting object can be set as the \code{response} argument
#' in \code{\link{dc_CA}} so as to give the same output as
#' a call to \code{\link{dc_CA}} with the abundance data as \code{response},
#' at least up to sign changes of the axes.
#'
#'

#' @inheritParams dc_CA
#' @param minimal_output logical. Default \code{TRUE} for use of the return value
#' as \code{response} in a call to \code{\link{dc_CA}}.
#'
#' @param verbose logical for printing a simple summary (default: TRUE)
#
#' @details
#' The argument \code{formulaTraits} determines which CWMs are calculated.
#' The CWMs are calculated for centered but further unstandardized trait data.
#' The argument \code{formulaEnv} determines which SNCs are calculated.
#' The SNCs are calculated for centered but further unstandardized environmental data.
#'
#' Empty (all zero) rows and columns in \code{response} are removed from the \code{response} and the corresponding
#' rows from \code{dataEnv} and \code{dataTraits}. Subsequently, any columns with missing values
#' are removed from  \code{dataEnv} and \code{dataTraits}. It gives an error (object 'name_of_variable' not found),
#' if variables with missing entries are specified in \code{formulaEnv} and \code{formulaTraits}.
#'
#' In the current implementation, \code{formulaEnv} and \code{formulaTraits} should
#' contain variable names as is, \emph{i.e.} transformations of variables in the formulas gives
#' an error ('undefined columns selected') when the \code{\link{scores}} function is applied.
#'
#' @returns
#'
#' The default returns a list of CWMs, SNCs, weights,
#' formulaTraits and a list of data with elements dataEnv and dataTraits.
#' When \code{minimal_output = FALSE}, many more statistics are given that are mainly technical or recomputed
#' when the return value is used as \code{response} in a call to \code{\link{dc_CA}}.
#'
#' @references
#' Kleyer, M., Dray, S., Bello, F., Lepš, J., Pakeman, R.J., Strauss, B., Thuiller,
#' W. & Lavorel, S. (2012) Assessing species and community functional responses to
#' environmental gradients: which multivariate methods?
#' Journal of Vegetation Science, 23, 805-821.
#' http://dx.doi.org/10.1111/j.1654-1103.2012.01402.x
#'
#' ter Braak, CJF, Šmilauer P, and Dray S. 2018. Algorithms and biplots for
#' double constrained correspondence analysis.
#' Environmental and Ecological Statistics, 25(2), 171-197.
#' https://doi.org/10.1007/s10651-017-0395-x or
#' http://rdcu.be/ETPh
#'
#' ter Braak C.J.F. and  P. Šmilauer  (2018). Canoco reference manual
#' and user's guide: software for ordination (version 5.1x).
#' Microcomputer Power, Ithaca, USA, 536 pp.
#'
#' Oksanen, J., et al. (2022)
#' vegan: Community Ecology Package. R package version 2.6-4.
#' http://CRAN.R-project.org/package=vegan.
#'
#' @seealso \code{\link{plot_dcCA}}, \code{\link{scores.dcca}}, \code{\link{print.dcca}} and \code{\link{anova.dcca}}
#' @example demo/dune_dcCA_via_CWMSNC.R
#' @export

fCWM_SNC <- function( response =NULL, dataEnv=NULL, dataTraits= NULL,
          formulaEnv = ~., formulaTraits = ~.,  divide.by.site.totals = TRUE, dc_CA_object  = NULL, minimal_output = TRUE, verbose = TRUE) {
  # response matrix or data frame, dataEnv and dataTraits data frames in which formualaE and formulaT are evaluated
  #dc_CA_object = result (value) of a previous run, can be used to save computing time for
  # runs that modify the formula for samples (step2: RDAonEnv) only
  # The step1 (CCAonTraits and the data and formulaTraits) are taken from dc_CA_object into the new result.
  # If set, formulaTraits, response, dataEnv, dataTraits are not used at all and have no efffect on the result
  call <- match.call()
 # if (is.null(dc_CA_object)){
    #  check and amend: make sure there are no empty rows or columns -----------------------------------------------------------------------
    if (any(is.na(response)))stop("The response should not have missing entries")
    if (any(response <0)) stop("The response should not have negative values")
    if (is.null(dataTraits)) stop("dataTraits must be specified in dc_CA")
    if (!is.matrix(response)) response <- as.matrix(response)

    id0 <-1
    while(length(id0)){
      TotR <- rowSums(response)
      id0 <- which(TotR == 0)
      if (length(id0)){
        response <- response[-id0,]
        dataEnv  <- dataEnv[-id0,]
      }
      TotC <- colSums(response)
      id0 <- which(TotC == 0)
      if (length(id0)){
        response <- response[,-id0]
        dataTraits <- dataTraits[-id0, ]
      }
    }

    # delete columns with missing data
    id = rep(FALSE, ncol(dataEnv))
    for (ii  in seq_along(id)){
      id[ii] <- sum(is.na(dataEnv[,ii]))==0
      if (!id[[ii]]) warning(
        paste("variable", names(dataEnv)[ii], "has missing values and is deleted from the environmental data")
      )
    }
    dataEnv <- dataEnv[, id]

    id = rep(FALSE, ncol(dataTraits))
    for (ii  in seq_along(id)){
      id[ii] <- sum(is.na(dataTraits[,ii]))==0
      if (!id[[ii]]) warning(
        paste("variable", names(dataTraits)[ii], "has missing values and is deleted from trait data")
      )
    }
    dataTraits <- dataTraits[, id]

    dataEnv <- as.data.frame(lapply(dataEnv, function(x){if (is.character(x)) x<- as.factor(x) else x; return(x) } ))
    dataTraits <- as.data.frame(lapply(dataTraits, function(x){if (is.character(x)) x<- as.factor(x) else x; return(x) } ))

    rownames(dataEnv) <- rownames(response)
    rownames(dataTraits) <- colnames(response)


    # end of check -----------------------------------------------------------------------
    if (divide.by.site.totals){
      response <- sweep(response, 1, STATS = TotR, FUN = '/')

    }
    Y <- as.matrix(response)/ sum(response)
    TotR <- rowSums(Y)
    TotC <- colSums(Y)
    weights = list(rows = TotR, columns = TotC) # unite sums
    Nobs = nrow(Y)
    # CWM and CWM ortho -------------------------------------------------------
   # formula = formulaTraits; data = dataTraits; w = weights$columns
    msqr <- msdvif(formulaTraits, dataTraits, weights$columns, XZ = FALSE)
    sWn <- sqrt(weights$columns)
    X = msqr$Xw/sWn; msd <- msqr$meansdvif
    traits2T_ortho <- qr.R(msqr$qrX)
    CWM2CWM_ortho <- solve(qr.R(msqr$qrX))
    # so Q represents the orthonormalized predictors
    CWM <- diag(1/weights$rows)%*% Y %*% X
    rownames(CWM) <- rownames(response)
    # # check
    #CWM_ortho <- CWM%*%CWM2CWM_ortho
    if (!minimal_output) T_ortho <- qr.Q(msqr$qrX)/ sWn
    #CWMs_orthonormal_traits <- diag(1/weights$rows)%*% Y %*% T_ortho
    #range(CWM_ortho-CWMs_orthonormal_traits)

    CWMs_orthonormal_traits <- CWM %*% CWM2CWM_ortho

    traits_explain <- sum(CWMs_orthonormal_traits^2 * weights$rows )
    #  else  traits_explain <- sum(CWMs_orthonormal_traits[,-1]^2 * weights$rows )

# SNC and SNC ortho -------------------------------------------------------
    #formula = formulaEnv; data = dataEnv; w = weights$rows
    msqr <- msdvif(formulaEnv, dataEnv, weights$rows, XZ = FALSE)
    sWn <- sqrt(weights$rows)
    X = msqr$Xw/sWn; msd <- msqr$meansdvif
    env2T_ortho <- qr.R(msqr$qrX)
    SNC2SNC_ortho <- solve(qr.R(msqr$qrX))
    # so Q represents the orthonormalized predictors
    SNC <- diag(1/weights$columns)%*% t(Y) %*% X
    rownames(SNC) <- colnames(response)
    SNCs_orthonormal_env <- SNC %*% SNC2SNC_ortho
    env_explain <- sum(SNCs_orthonormal_env^2 * weights$columns )
    # range(CWM1-CWMs_orthonormal_traits)# not OK!
# range(Xw -  qr.Q(qrX)%*% qr.R(qrX)) #0
# range(X -  Q%*% RfromQR) #
 if(minimal_output)
   out <- list(
     CWM = CWM,
     SNC = SNC,
     formulaEnv = formulaEnv,
     formulaTraits= formulaTraits, env_explain= env_explain,
     weights = weights, call = call, data = list(dataEnv = dataEnv, dataTraits = dataTraits)
   ) else out <- list(
      CWM = CWM, CWMs_orthonormal_traits = CWMs_orthonormal_traits,
      SNC = SNC, SNCs_orthonormal_env = SNCs_orthonormal_env,
      Nobs = Nobs, traits_explain = traits_explain, env_explain= env_explain,
      formulaEnv = formulaEnv, formulaTraits= formulaTraits,
      trans2ortho = list(CWM2CWM_ortho=CWM2CWM_ortho ,SNC2SNC_ortho=SNC2SNC_ortho),
      T_ortho = T_ortho,
      #E_ortho = E_ortho, T_ortho= T_ortho, Nobs = Nobs, divide.by.site.totals = divide.by.site.totals,RfromQR = RfromQR,
      weights = weights, call = call, data = list(dataEnv = dataEnv, dataTraits = dataTraits)
    )

  return(out)
}


f_canonical_coef_traits2 <- function(out){
  #B1_traitsN <- vegan::scores(step1,display = "reg", scaling = "sites", choices = 1:Rank_mod(step1))
  #B_star <- vegan::scores(step2, display = "species", scaling = "species", choices = 1:Rank_mod(step2))

  # do first the trivia of mean, sds and VIF
  # analogously  to calculate_b_se_tval in f_dc_CA.r
  if (is.null(out$CWM2CWM_ortho)) XZ= TRUE else XZ = FALSE
  ms <- msdvif(formula = out$formulaTraits, data=out$data$dataTraits, weights=out$weights$columns, XZ =XZ)
  avg = ms$meansdvif[,1]; sds = ms$meansdvif[,2]; VIF = ms$meansdvif[,3]
# end of trivia of mean sds and VIF
  if (!is.null(out$CWM2CWM_ortho))  {
    B1_traitsN <- out$CWM2CWM_ortho * sds
  }  else {
    if("cca" == class(out$CCAonTraits)[1]){
      B1_traitsN <- vegan::scores(out$CCAonTraits,display = "reg", scaling = "species", choices = 1:Rank_mod(out$CCAonTraits))
    } else {
      #  warning("Neither abundance data nor SNCs with transformation found; trait-related scores not available.")
      return(ms$meansdvif)
    }
  }

  if ("wrda" %in% class(out$RDAonEnv))  {
    B_star <- scores(out$RDAonEnv, display = "species", scaling = "sites", choices = 1: Rank_mod(out$RDAonEnv))
  } else  {
    B_star <- vegan::scores(out$RDAonEnv, display = "species", scaling = "sites", choices = 1:Rank_mod(out$RDAonEnv), const = 1)
  }

  c_traitsN <-  B1_traitsN  %*% B_star # End appendix 6.2 in ter Braak  et al 2018
  colnames(c_traitsN) <- paste("Regr", seq_len(ncol(c_traitsN)), sep= "")
  coef_normed <- cbind(Avg = avg , SDS = sds, VIF = VIF, c_traitsN, tval1 = NA)

  if (is.na(avg[1])){
    coef_normed[,"SDS"] <- NA
    attr(coef_normed, "warning")<-paste("regression coefs can only be interpreted columnwise in terms of their sign",
                                        "and quantitatively only rowwise", "as the standard deviation of traits is NA." )
  } else {
    attr(coef_normed, "meaning") <- "mean, sd, VIF, standardized regression coefficients, t-values not available"
  }
  return(coef_normed)
}

f2_orth <- function(CWM,formulaTraits, dataTraits, weights){

  msqr <- msdvif(formulaTraits, dataTraits, weights, XZ = FALSE)
  sWn <- sqrt(weights)
  X = msqr$Xw/sWn; msd <- msqr$meansdvif
  CWM2CWM_ortho <- solve(qr.R(msqr$qrX))
  colnames(CWM2CWM_ortho)<- rownames(CWM2CWM_ortho)
  CWMs_orthonormal_traits <- as.matrix(CWM)%*%CWM2CWM_ortho
  rownames(CWMs_orthonormal_traits) <- rownames(CWM)
  return(list(CWMs_orthonormal_traits=CWMs_orthonormal_traits, CWM2CWM_ortho=CWM2CWM_ortho))
}
