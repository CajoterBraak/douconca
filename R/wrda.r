#' @title Performs a weighted redundancy analysis
#'
#' @description
#' \code{wrda} is formula-based implementation of weighted redundancy analysis.
#'
#'
#' @param formula one or two-sided formula for the rows (samples) with row predictors in \code{data}.
#' The left hand side of the formula is ignored as it is specified in the next argument (\code{response}).
#' Specify row covariates (if any ) by adding \code{+ Condition(covariate-formula)}
#' to \code{formula} as in \code{\link[vegan]{rda}}. The \code{covariate-formula} should not contain a \code{~} (tilde).
#' @param response matrix or data frame of the abundance data (dimension \emph{n} x \emph{m}).
#' Rownames of \code{response}, if any, are carried through.
#' @param data matrix or data frame of the row predictors, with rows corresponding to those in \code{response}.
#' (dimension \emph{n} x \emph{p}).
#' @param weights row weights (a vector). If not specified unit weights are used.
#'
#' @param verbose logical for printing a simple summary (default: TRUE)
#
#' @details
#' In the current implementation, \code{formula} should
#' contain variable names as is, \emph{i.e.} transformations of variables in the formulas gives
#' an error ('undefined columns selected') when the \code{\link{scores}} function is applied.
#'
#' Compared to  \code{\link[vegan]{rda}}, \code{wrda} does not have residual axes,
#' (\emph{i.e.} no SVD or PCA of the residuals)
#'
#'
#'
#' @returns

#' All scores in the \code{dcca} object
#' are in scaling \code{"sites"} (1): the scaling with \emph{Focus on Case distances} .

#'
#' @references
#'
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
#' @example demo/dune_wrda.R
# @noRd
#' @export

wrda <- function(formula, response, data, weights= rep(1, nrow(data)), verbose = TRUE) {
  call <- match.call()
  Wn <- weights / sum(weights)
  sWn <- sqrt(Wn)
  msqr <- msdvif(formula, data=data, weights=Wn, XZ = TRUE)
  # transform to the unweighted case
  Yw <- as.matrix(response)*sWn
  # center
  Yw <- unweighted_lm_Orthnorm(Yw, matrix(sWn))
  total_variance <- sum(Yw^2)
  eY <- qr.resid(msqr$qrZ,Yw)
  Yfit_X <- qr.fitted(msqr$qrXZ,eY)
  svd_Yfit_X <- SVDfull(Yfit_X)
  # [1] "colsum"         "tot.chi"
  # [3] "Ybar"           "method"
  # [5] "call"           "pCCA"
  # [7] "CCA"            "CA"
  # [9] "inertia"        "regularization"
  # [11] "terms"          "terminfo"
  # needed?
  biplot <- NULL
  Nobs <- nrow(Yw)
  fct <- (Nobs/(Nobs-1))
  total_variance <-total_variance *fct
  ssY_gZ <- sum(eY^2) *fct
  ssY_XgZ <- sum(Yfit_X^2)*fct
  eig <- svd_Yfit_X$d^2 * (Nobs/ (Nobs -1))

  #step2$site_axes$site_scores$lc_env_scores <- step2$site_axes$site_scores$lc_env_scores * fct
  #step2$site_axes$site_scores$site_scores_unconstrained <- step2$site_axes$site_scores$site_scores_unconstrained * fct

  names(eig) <- paste("wRDA", seq_along(eig),sep="")
  CCA <- with(svd_Yfit_X, list(eig = eig, poseig = eig, u= u, v= v, rank= rank, qrank = msqr$qrXZ$rank,
                               tot.chi = ssY_XgZ, QR = msqr$qrXZ, biplot = biplot,
                               envcentre = NULL, centroids = NULL
                               ))
  if (ncol(msqr$Zw) ==1) pCCA <- NULL else
    pCCA  <- list(rank = min(ncol(Yw),msqr$qrZ$rank),tot.chi = total_variance - ssY_gZ,
                  QR = msqr$qrZ, envcentre =NULL)
  if (length(svd_Yfit_X$d)==1) diagd <-  matrix(svd_Yfit_X$d) else diagd<- diag(svd_Yfit_X$d)

  # need to be orthogonalized w.r.t Z
  #site_scores_unconstrained = qr.resid(qrZ, eY %*% svd_Yfit_X$v)/ sWn

  fct <- sqrt(fct)
  site_axes <- list(site_scores = list(
    site_scores_unconstrained =  qr.resid(msqr$qrZ, eY %*% svd_Yfit_X$v)/ (sWn/fct),
    lc_env_scores             =  (svd_Yfit_X$u %*% diag(svd_Yfit_X$d)) / (sWn/fct)
                                       )
  )

  species_axes <- list(species_scores = list(species_scores_unconstrained =svd_Yfit_X$v))

  object <- list(call = call, method = "wrda", tot.chi = total_variance,
           formula = formula, site_axes = site_axes, species_axes= species_axes,
           Nobs = Nobs, eigenvalues = eig,
           weights = list(rows = Wn, columns= rep(1/ncol(eY), ncol(eY)) ),
           data = data, eY = eY,
           pCCA = pCCA, CCA = CCA, CA = list(tot.chi = ssY_gZ-ssY_XgZ, rank = 0)

  )
  class(object) <- "wrda"
  #object$site_axes <- f_env_axes(object)
  #class(object) <- c("wrda",  "cca","dcca", class(object))

  return(object)
}

