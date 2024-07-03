#' @title Permutation Test for weighted redundancy analysis
#'
#' @description
#' \code{anova.wrda} performs residual predictor permutation (ter Braak 2022), 
#' which is robust against differences in the weights in weighted redundancy 
#' analysis (ter Braak, 2022). The arguments of the function are similar to 
#' those of \code{\link[vegan]{anova.cca}}, but more restricted.
#
#' @param object an object from \code{\link{dc_CA}}.
#' @param ... unused.
#' @param permutations a list of control values for the permutations as
#' returned by the function \code{\link[permute]{how}}, or the number of 
#' permutations required (default 999), or a permutation matrix where each row
#' gives the permuted indices.
#'
#' @param  by character \code{"axis"} which sets the test statistic to the
#' first eigenvalue of the dc-CA model. Default: \code{NULL} which set the test
#' statistic to the inertia named \code{constraintsTE} in the inertia element 
#' of \code{\link{dc_CA}}). This is the environmentally constrained inertia 
#' explained by the traits (without trait covariates). (which is equal to the 
#' trait-constrained inertia explained by the environmental predictors 
#' (without covariates).) The default is quicker computationally as it avoids
#' computation of an svd of permuted data sets.
#' 
#' @details
#' The algorithm is based on published R-code for residual predictor 
#' permutation in weighted redundancy analysis (ter Braak, 2022).
#'
#' @return
#' A list with two elements with names \code{table} and \code{eigenvalues}.
#' The \code{table} is as from \code{\link[vegan]{anova.cca}} and 
#' \code{eigenvalues} gives the dc-CA eigenvalues.
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
      stop("Error:: each row of permutations should have", N, "elements")
    }
  } else {
    stop("argument permutations should be integer, matrix or ", 
         "specified by permute::how().")
  }
  # Perform a weighted RDAR(M^*~E): an RDA of M^* on the environmental variables
  # using row weights R.
  sWn <- sqrt(object$weights$rows)
  Yw <-  object$eY * sqrt(object$Nobs/ (object$Nobs-1))
  msqr <- msdvif(object$formula, object$data, object$weights$rows, XZ = FALSE)
  Zw <- msqr$Zw 
  Xw <- msqr$Xw
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
  # what the env. variables explain of the trait-structured variation
  ss <- c(sapply(out_tes, function(x) x$ss[1]), out_tes[[length(out_tes)]]$ss[2])
  if (by == "axis") {
    df <- c(rep(1, length(ss) - 1), out_tes[[length(out_tes)]]$df[2])
    names(df) <- c(paste0("dcCA", seq_along(out_tes)), "Residual")
  } else {
    df <- out_tes[[length(out_tes)]]$df
    names(df) <- c("dcCA","Residual")
  }
  fraqExplained <- c(sapply(out_tes, function(x)x$ss[1]) / sum(out_tes[[1]]$ss), NA)
  F0 <- c(sapply(out_tes, function(x)x$F0[1]), NA)
  F.perm <- out_tes[[1]]$Fval
  R2.perm <- out_tes[[1]]$R2.perm
  if (length(out_tes)>1) {
    for (k in seq_along(out_tes)[-1]) {
      F.perm <- cbind(F.perm, out_tes[[k]]$Fval)
      R2.perm <- cbind(R2.perm,  out_tes[[k]]$R2.perm)
    }
  }
  p_val_axes1 <- c(cummax(sapply(out_tes, function(x)x$pval[1])), NA)
  eig <- out_tes[[1]]$eig
  names(eig)<- paste0("wRDA", seq_along(eig))
  axsig_dcCA_sites <- data.frame(df = df, Variance = ss, 
                                 R2 = fraqExplained, F = F0, 
                                 "Pr(>F)" = p_val_axes1, check.names = FALSE)
  object1 <- paste("Model:", c(object$call), "\n")
  heading <- paste0("Permutation test for weighted reduncancy analysis\n",
                    object1,
                    "Residualized predictor permutation\n",
                    howHead(attr(out_tes[[1]],"control") ))
  f_sites <- structure(axsig_dcCA_sites, heading = heading, 
                       control = attr(out_tes[[1]], "control"),
                       Random.seed = attr(out_tes[[1]], "seed"),
                       control = attr(out_tes[[1]], "control"),
                       F.perm = F.perm,
                       R2.perm = R2.perm,
                       class = c("anova.cca", "anova", "data.frame"))
  result <- list(table = f_sites, eig = eig)
  return(result)
}
