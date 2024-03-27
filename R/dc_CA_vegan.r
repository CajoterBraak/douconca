#' @title Performs a vegan-based double constrained correspondence Analysis (dc-CA)
#'
#' @description
#' \code{dcCA-vegan} performs double constained correspondence analysis (dc-CA) by
#' the two-step algorithm of ter Braak et al. (2018). This algorithm
#' combines and extends community- (sample-) and species-level analyses.
#' The first step uses \code{\link[vegan]{cca}} (Oksanen et al. 2022)
#' to regress the transposed abundance data on to the traits
#' and \code{\link[vegan]{rda}} to regress the community-weighted means (CWMs)
#' of the ortho-normalized traits on to the environmental predictors.
#' In this \code{vegan}-based version, the abundance data are divided by the sample total,
#' i.e. 'closed' in the terminology of compositional data analysis (CoDa). This
#' has the advantage that the results of a dc-CA correspond with unweighted (multi-trait)
#' community-level analyses, instead of corresponding to a weighted analysis (dc-CA weighs the
#' rows by the sample total, which is 1 after closure).
#' The current vegan-based analysis is efficient for community-level (site-based) permutation tests
#' but a factor 20 or so slower for species-level permutation tests.
#'
#' @param formulaEnv formula or one-sided formula for the rows (samples) with row predictors in \code{dataEnv}.
#' When two-sided, the left hand side of the formula is not used. Specify row covariates (if any ) by adding \code{+ Condition(covariate-formula)}
#' to \code{formulaEnv} as in \code{\link[vegan]{rda}}. Default: \code{~.}, i.e. all variables in \code{dataEnv} are predictor variables.
#' @param formulaTraits  formula or one-sided formula for the columns (species) with colum predictors in \code{dataTraits}.
#' When two-sided, the left hand side of the formula is not used. Specify column covariates (if any ) by  adding \code{+ Condition(covariate-formula)}
#' to \code{formulaTraits} as in \code{\link[vegan]{cca}}.Default: \code{~.}, i.e. all variables in \code{dataTraits} are predictor traits.
#' @param response matrix or data frame of the abundance data (dimension \emph{n} x \emph{m}).
#' Rownames of \code{response}, if any, are carried through.
#' @param dataEnv matrix or data frame of the row predictors, with rows corresponding to those in \code{response}.
#' ((dimension \emph{n} x \emph{p}).
#' @param dataTraits matrix or data frame of the column predictors,
#'  with rows corresponding to the columns in \code{response}.((dimension \emph{m} x \emph{q}).
#' @param dc_CA_vegan_object  optional object from an earlier run of this function. Useful if the
#' same formula for the columns (\code{formulaTraits}), \code{dataTraits} and \code{response} are used
#' with a new formula for the rows. If set, the data of the previous run is used and the result of its first step
#' is taken for the new analysis.
#' @param verbose logical for printing a simple summary (default: TRUE)
#
#' @details
#' Empty (all zero) rows and columns in \code{response} are removed from the \code{response} and the corresponding
#' rows from \code{dataEnv} and \code{dataTraits}. Subsequently, any columns with missing values
#' are removed from  \code{dataEnv} and \code{dataTraits}. It gives an error (object 'name_of_variable' not found),
#' if variables with missing entries are specified in \code{formulaEnv} and \code{formulaTraits}.
#' After 'closure' (division of the values in \code{response}
#' by their row sums), the subsequent algorithm follows the two-step algorithm of ter Braak et al. (2018).
#' and consits of two steps. First, the transpose of the \code{response} is regressed on to the traits
#' (the column predictors) using \code{\link[vegan]{cca}}
#' with \code{formulaTraits}.
#' The column scores of this analysis (in scaling 1) are community weigthed means (CWM) of the
#' orthonormalized traits.
#' These are then regressed on the environmental (row) predictors using \code{\link[vegan]{rda}} with
#' with \code{formulaEnv}.
#'
#'  All row based (sample-based) subsequent analyses, sample scores and permutation tests
#' can be obtained by the appropriate functions with argument \code{object$RDAonEnv}.
#'
#' The reason for the closure are
#' two-fold: dc-CA is the efficient summary of CWM-based regression analyses, which are rarely weighted, and
#' \code{vegan} \code{\link[vegan]{rda}} cannot do a weighted analysis, whereas \code{\link[vegan]{cca}} uses
#' the weights implied by the \code{response} after closure.
#'
#' The statistics and scores in the example \code{dune_dcCA.r},
#'  have been checked against the results in Canoco 5.15 (ter Braak & Smilauer, 1918).
#'
#' In the current implementation, \code{formulaEnv} and \code{formulaTraits} should
#' contain variable names as is, \emph{i.e.} transformations of variables in the formulas gives
#' an error ('undefined columns selected') when the \code{\link{scores}} function is applied.
#'
#' @returns
#' A list of \code{class} \code{dccav}; that is a list with elements
#' \describe{
#' \item{CCAonTraits}{a \code{\link[vegan]{cca.object}} from the \code{\link[vegan]{cca}} analysis
#' of the transpose of the closed \code{response} using formula \code{formulaTraits}. }
#' \item{formalaTraits}{the argument \code{formulaTraits}. }
#' \item{data}{a list of \code{Y} (response data after removing empty rows and columns and after closure)
#' and \code{dataEnv} and \code{dataTraits}.}
#' \item{weights}{a list of unit-sum weights of row and columns.
#' The names of the list are \code{c("row","columns")}, in that order.}
#' \item{Nobs}{number of sites (rows).}
#' \item{RDAonEnv}{a \code{\link[vegan]{cca.object}} from the \code{\link[vegan]{rda}} analysis
#' of the column scores of the \code{cca}, which are the CWMs of orthonormalized traits,
#' using formula \code{formulaEnv}. }
#' \item{formalaEnv}{the argument \code{formulaEnv}. }
#' \item{eigenvalues}{the dc-CA eigenvalues (same as those of the \code{\link[vegan]{rda}} analysis)}
#' \item{inertia}{a vector with four inertias (weighted variances):
#' \itemize{
#' \item total: the total inertia.
#' \item conditionT: the inertia explained by the condition in \code{formulaTraits} (neglecting row constraints).
#' \item traits_explain: the inertia explained by the traits (neglecting the row predictors and any
#' condition in \code{formulaTraits}).
#' This is the maximum that the row predictors could explain in dc-CA
#' (the sum of the following two items is thus less than this value).
#' \item conditionE: the trait-constrained inertia explained by the condition in \code{formulaEnv}.
#' \item constraintsTE: the trait-constrained inertia explained by the predictors (without the row covariates).
#' }
#'  }
#' }
#' If \code{verbose} is \code{TRUE} (or after \code{out <-print(out)} is invoked )
#' there are three more items (in this version).
#' \itemize{
#' \item \code{c_traits_normed}: mean, sd, VIF and (regression) coefficients of
#'  the traits that define the dc-CA trait axes (composite traits), and their optimistic t-ratio.
#' \item \code{c_env_normed}:  mean, sd, VIF and (regression) coefficients of the environmental variables that define the dc-CA axes
#'  in terms of the environmental variables (composite gradients), and their optimistic t-ratio.
#' \item \code{species_axes}: a list with four items
#'  \itemize{
#'  \item \code{species_scores}: a list with names \code{c("species_scores_unconstrained", "lc_traits_scores")} with the
#'  matrix with species niche centroids along the dc-CA axes (composite gradients) and
#'  the matrix with linear combinations of traits.
#'  \item \code{correlation}: a matrix with inter-set correlations of the traits with their SNCs.
#'  \item \code{b_se}: a matrix with (unstandardized) regression coefficients for traits and their optimistic standard errors.
#'  \item \code{R2_traits}: a vector with coefficient of determination (R2) of the SNCs on to the traits.
#'  The square-root thereof could be called the species-trait correlation in analogy with
#'  the species-environment correlation in CCA.
#'  }
#'  \item \code{sites_axes}: a list with four items
#'  \itemize{
#'  \item \code{site_scores}: a list with names \code{c("site_scores_unconstrained", "lc_env_scores")} with the
#'  matrix with community weighted means (CWMs) along the dc-CA axes (composite gradients) and
#'  the matrix with linear combinations of environmental variables.
#'  \item \code{correlation}: a matrix with inter-set correlations of the environmental variables with their CWMs.
#'  \item \code{b_se}: a matrix with (unstandardized) regression coefficients for environmental
#'  variables and their optimistic standard errors.
#'  \item \code{R2_env}: a vector with coefficient of determination (R2) of the CWMs on to the environmental variables.
#'  The square-root thereof has been called the species-environmental correlation in CCA.
#'  }
#'
#' }
#' All scores in the \code{dccav} object
#' are in scaling \code{"sites"} (1): the scaling with \emph{Focus on Case distances} .

#'
#' @references
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
#' @seealso \code{\link{plot_dcCA}}, \code{\link{scores.dccav}}, \code{\link{print.dccav}} and \code{\link{anova.dccav}}
#' @example demo/dune_dcCA.R
#' @export

dc_CA_vegan <- function(formulaEnv = ~., formulaTraits = ~., response =NULL, dataEnv, dataTraits= NULL, dc_CA_vegan_object  = NULL, verbose = TRUE) {
  # response matrix or data frame, dataEnv and dataTraits data frames in which formualaE and formulaT are evaluated
  #dc_CA_vegan_object = result (value) of a previous run, can be used to save computing time for
  # runs that modify the formula for samples (step2: RDAonEnv) only
  # The step1 (CCAonTraits and the data and formulaTraits) are taken from dc_CA_vegan_object into the new result.
  # If set, formulaTraits, response, dataEnv, dataTraits are not used at all and have no efffect on the result

  if (is.null(dc_CA_vegan_object)){
    #  check and amend: make sure there are no empty rows or columns -----------------------------------------------------------------------
    if (any(is.na(response)))stop("The response should not have missing entries")
    if (any(response <0)) stop("The response should not have negative values")
    if (is.null(dataTraits)) stop("dataTraits must be specified in dc_CA_vegan")
    if (!is.matrix(response)) response <- as.matrix(response) else stop("response (matrix or df) must specified")
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
    # close the data (divide by the row total, to get strictly compositional data) -----------------------------------------------------------------------

    call <- match.call()
    Abun_frac <- sweep(response, 1, STATS = TotR, FUN = '/')
    TotRfrac <- rowSums(Abun_frac)
    TotC <- colSums(Abun_frac)
    tY <- t(Abun_frac)
    formulaTraits <- change_reponse(formulaTraits, "tY")
    environment(formulaTraits)<- environment()
    step1 <-vegan::cca(formulaTraits, data = dataTraits)
    data= list(Y = Abun_frac, dataEnv = dataEnv, dataTraits = dataTraits)
    out1 <- list(CCAonTraits = step1,
                 formulaTraits= formulaTraits,
                 data = list(Y = Abun_frac, dataEnv = dataEnv, dataTraits = dataTraits),
                 call = call,
                 weights = list(columns = TotC/sum(TotC))
    )
  } else {
    step1 <- dc_CA_vegan_object$CCAonTraits
    out1 <- dc_CA_vegan_object[c("CCAonTraits", "formulaTraits","data", "weights","call")]
  }

  n <- nrow(out1$data$Y)
  CWMs_orthonormal_traits <- vegan::scores(step1, display= "species",
                scaling = "species",
                choices = seq_len(Rank_mod(out1$CCAonTraits)) ) * sqrt((n-1)/n)
  if (rownames(CWMs_orthonormal_traits)[1]=="col1") rownames(CWMs_orthonormal_traits) <- paste("Sam", seq_len((nrow(out1$data$dataEnv))),sep="")
  formulaEnv <- change_reponse(formulaEnv, "CWMs_orthonormal_traits")
  environment(formulaEnv)<- environment()
  step2 <- vegan::rda(formulaEnv, data = out1$data$dataEnv)



  out <- c(out1[-5], list(RDAonEnv = step2,
                      formulaEnv = formulaEnv,
                      eigenvalues =  vegan::eigenvals(step2, model = "constrained"),
                      weights = list(rows = rep(1/n, n),columns = out1$weights$columns),
                      Nobs = n
                      )
  )


  inertia <- try(f_inertia(out))

  if("try-error"%in% class(inertia)) {warning("could not obtain inertia's"); print(inertia)}

  out$inertia <- inertia


  if (verbose) {
    out<-print.dccav(out)
  }

  class(out) <- c("dccav", "list")


  return(out)
}

