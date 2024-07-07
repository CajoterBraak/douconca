#' @title Get constraining and conditional formulas from an rda or cca ordination
#'
#' @description
#' \code{get_Z_X_XZ_formula} derives the focal constraining formula and 
#' (if present) the conditioning formula from a vegan formula 
#' \code{~X+Condition(Z)}.
#' 
#' @param formula  a vegan formula of the form \code{~X+Condition(Z)} with
#' \code{X} and \code{Z} formulas (without a tilde). 
#' Dimension reduction is applied to the effects
#' of \code{X}  given the effects of \code{Z} on the response.
#' 
#' @param data matrix or data frame of the row predictors. 
#' Default \code{NULL} in which case the result values for elements
#'  "focal_factor" "Condi_factor" are '\code{NA}.
#' 
#' @return A list with elements
#'[1] "formula_XZ"   "formula_X0"   "formula_Z"    "focal_nams"  
#'[5] "Condi_nams"   "formula_X1"   "focal_factor" "Condi_factor"
#'[9] "all_nams" 
#' 
#' @details
#' In the current implementation, \code{formula} should contain variable names
#' as is, \emph{i.e.} transformations of variables in the formulas gives
#' an error ('undefined columns selected') when the \code{\link{scores}} 
#' function is applied. Interactions and nesting ("*", ":" and "/" ) are allowed.
#' 
#' #' 
#' @noRd
#' @keywords internal
get_Z_X_XZ_formula <- function(formula, 
                               data = NULL) {
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
