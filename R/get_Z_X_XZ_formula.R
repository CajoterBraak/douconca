#' @title Get constraining and conditional formulas from an rda or cca ordination
#'
#' @description
#' \code{get_Z_X_XZ_formula} derives the focal constraining formula and 
#' (if present) the conditioning formula from a vegan formula
#' \code{~X+Condition(Z)}.
#'  
#' @param formula  a vegan formula of the form \code{~X+Condition(Z)} with
#' \code{X} and \code{Z} formulas (without a tilde). Dimension reduction is 
#' applied to the effects of \code{X}  given the effects of \code{Z} on the 
#' response.
#' 
#' @param data matrix or data frame of the row predictors. Default \code{NULL},
#' in which case the result values for elements "focal_factor" "Condi_factor" 
#' are '\code{NA}.
#' 
#' @return A list with elements:
#' "formula_XZ", "formula_X0", "formula_Z", "focal_nams", "Condi_nams",
#' "formula_X1", "focal_factor", "Condi_factor", and "all_nams".
#' 
#' @details
#' In the current implementation, \code{formula} should contain variable names
#' as is, \emph{i.e.} transformations of variables in the formulas gives
#' an error ('undefined columns selected') when the \code{\link{scores}} 
#' function is applied. Interactions and nesting ("*", ":" and "/" ) are 
#' allowed.
#' 
#' @noRd
#' @keywords internal
get_Z_X_XZ_formula <- function(formula, 
                               data = NULL) {
  # get formulas and focal and conditioning_factors from a formula
  # with optional data to determine which variables are factors
  trms <- delete.response(terms(formula, specials = "Condition", 
                                keep.order = TRUE))
  trmLabs <- rownames(attr(trms, "factors"))
  condId <- attr(trms, "specials")$Condition
  if (!length(condId)) {
    condTrms <- NULL
    formula_Z <- ~1
  } else {
    condi <- trmLabs[condId]
    formula_Z <- as.formula(paste("~", substr(condi, 11, nchar(condi) - 1)))
    condTrms <- attr(terms(formula_Z), "term.labels")
    trmLabs <- trmLabs[-condId]
  }
  formula_X0 <- reformulate(trmLabs, intercept = FALSE)
  formula_X1 <- reformulate(trmLabs, intercept = TRUE)
  formula_XZ <- reformulate(c(condTrms, trmLabs), intercept = TRUE)
  focal_nams <- all.vars(formula_X0)
  condi_nams <- all.vars(formula_Z) 
  if (!is.null(data)) {
    mff <- model.frame(formula_X0, data = data)
    trmsf <- attr(mff, "terms")
    focal_factor <- all.vars(trmsf)[attr(trmsf, "dataClasses") == "factor"]
    mfc <- model.frame(formula_Z, data = data)
    trmsc <- attr(mfc, "terms")
    condi_factor <- all.vars(trmsc)[attr(trmsc, "dataClasses") == "factor"]
  } else {
    focal_factor <- condi_factor <- NA
  }
  return(list(formula_XZ = formula_XZ,
             # formula_X0 = formula_X0,
              formula_Z = formula_Z,
              focal_nams = focal_nams,
              Condi_nams = condi_nams,
              formula_X1 = formula_X1,
              focal_factor = focal_factor,
              Condi_factor = condi_factor,
              all_nams = union(focal_nams, condi_nams)
  ))
}
