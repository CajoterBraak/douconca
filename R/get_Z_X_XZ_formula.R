#' @title Get constraining and conditional formulas from an rda or cca ordination
#'
#' @description
#' \code{get_Z_X_XZ_formula} derives the focal constraining formula and 
#' (if present) the conditioning formula from a vegan formula.
#' 
#' @param objecta result of \code{\link[vegan]{rda}}, \code{\link[vegan]{cca}} 
#' specified via a formula (S3 method for class 'formula'), or a result of
#' \code{\link{doPRC}} or \code{\link{PRC_scores}}.
#' 
#' @return A list with element \code{`focal factor`} and \code{condition}
#' 
#' @example demo/PRC_pyrifos_bk.R
#' 
#' @references
#' ter Braak C.J.F. & Šmilauer P. (2018): Canoco reference manual and user's guide:
#' software for ordination, version 5.1x. Microcomputer Power, Ithaca, USA, 536 pp.
#' (\url{http::www.canoco5.com})
#' 
#' @seealso \code{\link{doPRC}}, \code{\link{PRC_scores}}
#' 
#' @noRd
#' @keywords internal
get_Z_X_XZ_formula <- function(formula, 
                               data = NULL, 
                               factors_only = FALSE) {
  # get formulas and focal and conditioning_factors from a formula
  # with optional data to determine which variables are factors
  tl <- attr(terms(formula, data = data), "term.labels")
  if (is.null(tl)) {
    out <- list(formula_XZ = NULL, formula_X = NULL, formula_Z = NULL)
    names(out) <- c("focal factor", "condition")
    stop("Specify wrda or cca via a formula, so that predictor and covariate ", 
         "formulas can be found.")
  }
  # find all variable names by deleting from tl Condition and all interactions
  # first delete Condition
  idC <- pmatch("Condition(", tl)
  if (!is.na(idC)) {
    # get formula of the Condition
    condi <- tl[idC]
    fC_wo_tilde <- substr(condi, 11, nchar(condi) - 1)
    fC_formula <- as.formula(paste("~", fC_wo_tilde))
  } else { 
    fC_formula <- ~1
  }
  if (!is.na(idC)) {
    tl <- tl[-idC]
  }
  fFocalX <- paste(tl, collapse = "+")
  if (!is.na(idC)) {
    fFocalXZ <- as.formula(paste("~", fC_wo_tilde, "+", fFocalX)) 
  } else {
    fFocalXZ <- as.formula(paste("~", fFocalX))
  }
  fFocalX0 <- as.formula(paste("~0 +", fFocalX))
  fFocalX1 <- as.formula(paste("~1 +", fFocalX))
  focal_nams <- unique(unlist(strsplit(tl, split = ":", fixed = TRUE)))
  Condi_nams <- unique(unlist(strsplit(attr(stats::terms(fC_formula), "term.labels"), 
                                       split = ":", fixed = TRUE)))
  focal_nams <- focal_nams[!focal_nams %in% Condi_nams]
  if (!is.null(data)) {
    fctrs <- names(data)[sapply(data, is.factor)]
    focal_factor <- focal_nams[focal_nams %in% fctrs]
    Condi_factor <- Condi_nams[Condi_nams %in% fctrs]
    if (length(focal_factor) == 0) {
      focal_factor <- NULL
    } 
    if (length(Condi_factor) == 0) {
      Condi_factor <- NULL
    }
  } else {
    focal_factor <- Condi_factor <- NA
  }
  
  return(list(formula_XZ = fFocalXZ, 
              formula_X0 = fFocalX0, 
              formula_Z = fC_formula,
              focal_nams = focal_nams,
              Condi_nams = Condi_nams, 
              formula_X1 = fFocalX1,
              focal_factor = focal_factor,
              Condi_factor = Condi_factor,
              all_nams = unique(c(focal_nams, Condi_nams))
  ))
}
