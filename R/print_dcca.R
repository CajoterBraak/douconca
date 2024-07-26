#' Print a summary of a dc-CA object
#' 
#' @param x a dc-CA object from \code{\link{dc_CA}}
#' @param ...  Other arguments passed to the function (currently ignored).
#' 
#' @details
#' \code{x <- print(x)} is more efficient for \code{\link{scores.dcca}} than 
#' just \code{print(x)} if \code{\link{dc_CA}} is called without argument
#' \code{verbose} (or called with \code{verbose = FALSE}).
#' 
#' @example demo/dune_dcCA.r
#' 
#' @noRd
#' @keywords internal
print_dcca <- function(x, 
                       ...) {
  if (inherits(x, "dcca")) {
    if (!"species_axes" %in% names(x)) {
      x$site_axes <- f_env_axes(x)
      x$species_axes <- f_trait_axes(x)
      x$c_env_normed <- x$site_axes$c_env_normed
      x$c_traits_normed <- x$species_axes$c_traits_normed
    }
  } else if (inherits(x, "wrda", which = TRUE) == 1) {
    x$site_axes <- f_env_axes(x)
    x$c_env_normed <- x$site_axes$c_env_normed
  } else {
    stop("The argument must be of class 'dcca' or 'wrda', the result of ", 
         "the function dc_CA or wrda.\n")
  }
  choices <- 1:4
  if (inherits(x, "dcca")) {
    cat("Step 1: the CCA ordination of the transposed matrix with", 
        "trait constraints,\n")
    cat("        useful in itself and also yielding CWMs of the", 
        "orthonormalized traits for step 2.\n")
    print(x$CCAonTraits)
    cat("Step 2: the RDA ordination of CWMs of the orthonormalized traits \n", 
        "       of step 1 with environmental constraints:\n")
    print(x$RDAonEnv)
    c_t <- x$c_traits_normed[, c(choices, 4 + rank_mod(x)), drop = FALSE]
  } else {
    print.cca(x)
    c_t <- NULL
  }
  if (any(inherits(x, c("dccav", "wrda"), which = TRUE) == 1)) {
    c_e <- x$c_env_normed[, c(choices, 4 + rank_mod(x)), drop = FALSE]
    cat("mean, sd, VIF and canonical coefficients with their optimistic [!]", 
        "t-values:\n")
    print(round(c_e, 4))
  }
  if (!is.null(c_t)) {
    print (round(c_t, 4))
  }
  cat("\n")
  if (inherits(x, "dcca")) {
    print(round(x$inertia, 3))
  }
  return(invisible(x))
}
