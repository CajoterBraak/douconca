#' Print a summary of a dc-CA object
#' @param x a dc-CA object from \code{\link{dc_CA_vegan}}
#' @param ...  Other arguments passed to the function (currently ignored).
#' @example demo/dune_dcCA.r
#' @details
#' \code{x <- print(x)} is more efficient for \code{\link{scores.dccav}} than just \code{print(x)}
#' if \code{\link{dc_CA_vegan}} is called without argument \code{verbose} (or called with \code{verbose = FALSE}).
#'
#' @export
#'
print.dccav <- function(x, ...){

  if (!"species_axes"%in%names(x)){
    x$site_axes <-  f_env_axes(x)
    x$species_axes <- f_trait_axes(x)
    x$c_env_normed <- x$site_axes$c_env_normed
    x$c_traits_normed <- x$species_axes$c_traits_normed
  }
  choices <- c(1:4);
  cat("Step 1: the CCA ordination of the transposed matrix with trait constraints,\n")
  cat("        useful in itself and also yielding CWMs of the orthonormalized traits for step 2.\n")

  print(x$CCAonTraits)
  cat("Step 2: the RDA ordination of CWMs of the orthonormalized traits \n        of step 1 with environmental constraints:\n")
  print(x$RDAonEnv)
  c_e <- x$c_env_normed[,c(choices, 4 + Rank_mod(x$RDAonEnv))]
  c_t <- x$c_traits_normed[,c(choices, 4 + Rank_mod(x$RDAonEnv))]

  cat("mean, sd, VIF and canonical coefficients with their optimistic [!] t-values:\n")
  print (round(c_e,4))

  print (round(c_t,4))
  #cat("\nWarning: the t-values are optimistic, i.e. an underestimate of their true absolute value")
  cat("\n")
  #inertia <- x$inertia[, 1:ncol(x$inertia),drop= FALSE]
  print (round(x$inertia,3))
  # expla <- c("total inertia","inertia of the trait condition", "trait-constrained inertia",
  #            "trait-constrained inertia explained by the condition in formulaEnv",
  #            "trait-constrained inertia explained by the predictors in formulaEnv")
  # names(expla) <- c("total","conditionT","traits_explain","conditionE","constraintsE")
  # print(matrix( expla[rownames(x$inertia)], ncol=1,
  #                                              dimnames = list(rownames(x$inertia),"meaning")))


  class(x) <- c("dccav", "list")
  invisible(x)
}
