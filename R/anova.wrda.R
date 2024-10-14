#' @title Permutation Test for weighted redundancy analysis
#'
#' @description
#' \code{anova.wrda} performs residual predictor permutation for weighted 
#' redundancy analysis (wRDA), which is robust against differences in the 
#' weights (ter Braak, 2022). The arguments of the function are similar to 
#' those of \code{\link[vegan]{anova.cca}}, but more restricted.
#
#' @param object an object from \code{\link{dc_CA}}.
#' @param ... unused.
#' @param permutations a list of control values for the permutations as
#' returned by the function \code{\link[permute]{how}}, or the number of 
#' permutations required (default 999), or a permutation matrix where each row
#' gives the permuted indices.
#'
#' @param by character \code{"axis"} which sets the test statistic to the
#' first eigenvalue of the RDA model. Default: \code{NULL} which sets the test
#' statistic to the weighted variance fitted by the predictors
#' (=sum of all constrained eigenvalues). The default is quicker 
#' computationally as it avoids computation of an svd of permuted data sets.
#' 
#' @details
#' The algorithm is based on published R-code for residual predictor 
#' permutation in weighted redundancy analysis (ter Braak, 2022), but
#' using QR-decomposition instead of ad-hoc least-squares functions.
#'
#' @return
#' A list with two elements with names \code{table} and \code{eigenvalues}.
#' The \code{table} is as from \code{\link[vegan]{anova.cca}} and 
#' \code{eigenvalues} gives the wrda eigenvalues.
#' 
#' @references
#' ter Braak, C.J.F. (2022) Predictor versus response permutation
#' for significance testing in weighted regression and redundancy analysis.
#' Journal of statistical computation and simulation, 92, 2041-2059.
#' \doi{10.1080/00949655.2021.2019256}
#' 
#' @example demo/dune_wrda.R
#' 
#' @importFrom stats anova
#' @export
anova.wrda <- function(object, 
                       ...,
                       permutations = 999, 
                       by = c("omnibus", "axis")) {
  # permat  a matrix of permutations. 
  # If set overrules permuations.
  by <- match.arg(by)
  N <- nrow(object$data) 
  if (inherits(permutations, c("numeric", "how", "matrix"))) {
    if (is.numeric(permutations) && !is.matrix(permutations)) {
      permutations <- permute::how(nperm = permutations[1])
    } else if (is.matrix(permutations) && ncol(permutations) != N) {
      stop("each row of permutations should have", N, "elements.\n")
    }
  } else {
    stop("argument permutations should be integer, matrix or ", 
         "specified by permute::how().\n")
  }
  # Perform a weighted RDAR(M^*~E): an RDA of M^* on the environmental variables
  # using row weights R.
  sWn <- sqrt(object$weights$rows)
  Yw <-  object$eY * sqrt(object$Nobs/ (object$Nobs-1))
  msqr <- msdvif(object$formula, object$data, object$weights$rows, XZ = FALSE)
  Zw <- msqr$Zw 
  Xw <- msqr$Xw
  dfpartial = msqr$qrZ$rank
  # residual predictor permutation.
  out_tes <- list()
  out_tes[[1]]  <- randperm_eX0sqrtw(Yw,Xw, Zw, sWn = sWn, 
                                     permutations = permutations, 
                                     by = by, return = "all")
  if (by == "axis") {
    while (out_tes[[1]]$rank > length(out_tes)) {
      Zw <- cbind(Zw, out_tes[[length(out_tes)]]$EigVector1)
      out_tes[[length(out_tes) + 1]] <- 
        randperm_eX0sqrtw(Yw,Xw, Zw, sWn = sWn, 
                          permutations = permutations, by = by, return = "all")
    }
  }
  f_sites <- fanovatable(out_tes, Nobs = N, dfpartial= dfpartial, type= "wrda", 
                         calltext= c(object$call))  
  result <- list(table = f_sites, eigenvalues = attr(f_sites, "eig"))
  return(result)
}
