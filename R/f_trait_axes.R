#' Utility function: obtain column (species) related Scores
#' 
#' @param out object from \code{\link{dc_CA}}
#' @param which_cor character or names of traits for which inter-set 
#' correlations must calculated. Default: "in_model" for all traits and 
#' variables in the model specified by \code{formulaTraits}.
#' 
#' @noRd
#' @keywords internal
f_trait_axes <- function(out, 
                         which_cor = "in model") {
  # which_cor character or names of traits
  # for which inter-set correlations must calculated.
  #  Default "in_model" for all traits and variables in the model
  # SNC lc_traits and trait regr, tval, cor
  if (is.null(out$data$Y) && is.null(out$SNC)) {
    warning("SNC analysis is not available (e.g. no unconstrained species scores).")
  }
  if (is.null(out$data$Y)) {
    # NB: need to check whether the axes have not been flipped
    #diag(wcor (mod$species_axes$species_scores$lc_traits_scores, SNC, w =  mod$weights$columns))
    # todo:not fully correct yet note that out$CCAonTraits may not exist and
    # out$RDAonEnv may be a wrda
    # NB: we need the rotation of the dataTraits to the orthonormalized traits of which
    # the CWMs are taken.
    if (!is.null(out$SNCs_orthonormal_env)) {
      SNC <- out$SNCs_orthonormal_env
      step2 <- wrda(formula = out$formulaTraits, 
                    response = SNC * sqrt((nrow(SNC) - 1) / nrow(SNC)),
                    data = out$data$dataTraits, weights = out$weights$columns)
      if (diff(range((step2$eigenvalues + 1.e-10) / 
                     (out$eigenvalues+1.e-10))) > 1.6e-3) {
        warning("\nThe eigenvalues of the CWM and SNC analyses differ.",
                "\nCWM analysis has eigenvalues\n", round(out$eigenvalues, 6),
                "\nSNC analysis has eigenvalues\n", round(step2$eigenvalues, 6))
      }
      SNC <- step2$site_axes$site_scores$site_scores_unconstrained
      out$CCAonTraits$CCA$QR <- step2$CCA$QR
    } else { #is.null(out$SNCs_orthonormal_env)
      if (!is.null(out$CWM2CWM_ortho)) {
        warning("Trait regression coefficients are derived from the CWM analysis")
        regr <- out$c_traits_normed0[, -c(1, 2, 3)]
        regr <- regr[, -ncol(regr), drop = FALSE]
        # must use model.matrix
        fX <- get_Z_X_XZ_formula(out$formulaTraits, out$data$dataTraits)
        X <- model.matrix(fX$formula_X1, 
                          data = out$data$dataTraits)[, -1, drop = FALSE]
        lc_trait_scores <- standardize_w(X) %*% regr
      } else  {
        warning("Trait regression coefficients  are not availabe")
        lc_trait_scores <-NULL
      }
      return(
        list(species_scores = list(species_scores_unconstrained = NULL,
                                   lc_traits_scores = lc_trait_scores), 
             correlation = NULL, 
             c_traits_normed = out$c_traits_normed0, 
             b_se = NULL, 
             R2_traits = NULL))
    }
  } else { # try to obtain SNCs
    step2 <- NULL
    if (!is.null(out$data$Y)) {
      if (inherits(out, "dccav")) {
        lc_scores  <- scores(out$RDAonEnv, display = "lc", 
                             scaling = "species",
                             choices = seq_len(Rank_mod(out$RDAonEnv)), 
                             const = sqrt(out$Nobs))
      } else if (inherits(out, "dcca")) {
        lc_scores <- out$RDAonEnv$site_axes$site_scores$lc_env_scores
      } else {
        stop("first argument must be of class dcca or dccav")
      }
      SNC <-  (t(as.matrix(out$data$Y)) %*% lc_scores) / 
        (out$weights$columns * out$Nobs)
    } else {
      warning("something wrong in f_trait_axes")
      return(NULL)}
  }
  if (!is.null(out$CCAonTraits$pCCA)) {
    # orthogalize with respect to any covariate
    SNC <- SNC - calculate_b_se_tval(get_QR(out$CCAonTraits, model = "pCCA"), 
                                     y = SNC, w = out$weights$columns,  
                                     scale2 = 0, name = "SNC", 
                                     fitted_only = TRUE)
  }
  res <- calculate_b_se_tval(get_QR(out$CCAonTraits), y = SNC, 
                             w = out$weights$columns,  scale2 = 1, name = "SNC")
  c_traits_normed <- res$coef_normed
  attr(c_traits_normed, which = "warning") <-
    "The t-values are optimistic, i.e. an underestimate of their true absolute value"
  # check sign of axes when !is.null(out$SNCs_orthonormal_env)
  if (!is.null(out$SNCs_orthonormal_env) && is.null(out$data$Y)) {
    if(!is.null(out$c_traits_normed0)) {
      rseq <- seq_len(Rank_mod(out))
      ncov <- nrow(c_traits_normed) - nrow(out$c_traits_normed0)
      if (ncov) {
        ratio <- sign(colSums(sign(c_traits_normed[-seq_len(ncov), 
                                                   3 + rseq, drop = FALSE] / 
                                     out$c_traits_normed0[,3 + rseq]))) 
      } else {
        ratio <- sign(colSums(sign(c_traits_normed[, 3 + rseq, drop = FALSE] / 
                                     out$c_traits_normed0[, 3 + rseq])))
      }
      if (length(ratio) == 1) {
        flip <- matrix(ratio) 
      } else {
        flip <- diag(ratio)
      }
      c_traits_normed[, 3 + rseq] <- 
        c_traits_normed[, 3 + rseq, drop = FALSE] %*% flip
      c_traits_normed[, 3 + length(rseq) + rseq] <-
        c_traits_normed[, 3 + length(rseq) + rseq, drop = FALSE] %*% flip
      res$y <- res$y %*% flip
      res$fitted <- res$fitted %*% flip
      res$b_se[, rseq] <- as.matrix(res$b_se[, rseq, drop = FALSE]) %*% flip
      SNC <- SNC %*% flip
    } else { 
      warning("The orientation of the trait and environmental axes cannot be ", 
      "aligned, as information is missing how the orthonormized traits must ", 
      "be backtransformed to the original traits.")
    }
  }
  # correlations of the dataTraits with the SNC wrt the axes
  if (which_cor[1] == "in model") {
    whichc <- get_Z_X_XZ_formula(out$formulaTraits, out$data$dataTraits)$focal_nams
    gg <- get_Z_X_XZ_formula(out$formulaTraits, out$data$dataTraits)
    traits0 <- model.matrix(gg$formula_X0, data = out$data$dataTraits)
  } else {
    whichc <- which_cor
    traits0 <- model.matrix(~.-1, 
                            data = out$data$dataTraits[, whichc, drop = FALSE])
  }
  Cor_Trait_SNC <- wcor(traits0, SNC, w = out$weights$columns)
  colnames(Cor_Trait_SNC) <- paste0("SNC-ax", seq_len(ncol(Cor_Trait_SNC)))
  attr(Cor_Trait_SNC, which = "meaning") <- 
    "inter set correlation, correlation between traits and SNC of axes"
  out2 <- list(species_scores = list(species_scores_unconstrained = res$y,
                                     lc_traits_scores = res$fitted), 
               correlation = Cor_Trait_SNC, 
               c_traits_normed = c_traits_normed, 
               b_se = res$b_se, 
               R2_traits = res$R2, 
               step2 = step2)
  return(out2)
}
