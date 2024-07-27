# utility functions for f_trait_axes and f_env_axes in print.dcca and scores.dcca--------------

#' @noRd
#' @keywords internal
wcor <- function(X, 
                 Y = X,
                 w = rep(1, nrow(X))) {
  # weighted correlation between matrix X and Y
  w <- w / sum(w)
  Xstd <- standardize_w(X, w)
  Ystd <- standardize_w(Y, w)
  t(Xstd) %*% diag(w) %*% Ystd
}

#' @noRd
#' @keywords internal
mean_w <- function(X,
                   w = rep(1 / nrow(X), nrow(X))) {
  t(w / sum(w)) %*% X
}

#' @noRd
#' @keywords internal
center_w <- function(X,
                     w = rep(1 / nrow(X), nrow(X))) {
  X - rep(1, length(w)) %*% t(w) %*% X 
}

#' @noRd
#' @keywords internal
standardize_w <- function(X, 
                          w = rep(1 / nrow(X), nrow(X))) {
  # NB requires w to be have sum 1
  ones <- rep(1, length(w))
  Xc <- X - ones %*% t(w) %*% X
  Xstd <- Xc / ones %*% sqrt(t(ones) %*% (Xc * Xc * w))
  return(Xstd)
}

#' @noRd
#' @keywords internal
mean_sd_w <- function(X,
                      w = rep(1 / nrow(X), nrow(X))){
  # NB requires w to be have sum 1
  ones <- rep(1, length(w))
  mean_w <- t(w / sum(w)) %*% X
  Xc <- X - ones %*% mean_w
  sd_w <-  sqrt(t(ones) %*% (Xc * Xc * w))
  return(list(mean = mean_w, sd = sd_w))
}

#' @noRd
#' @keywords internal
msdvif <- function(formula = NULL, 
                   data, 
                   weights, 
                   XZ = FALSE, 
                   novif = FALSE, 
                   contrast = TRUE) {
  # calc mean variance and vif from for X given Z or XZ 
  # with qr of X|Z or of centered XZ
  if (is.null(formula)) {
    f <- ~. 
  } else {
    f <- formula
  }
  ff <- get_Z_X_XZ_formula(f, data)
  if (XZ) {
    X <- model.matrix(ff$formula_XZ, data = data)[, -1, drop = FALSE]
  } else {
    if (contrast) {
      X <- model.matrix(ff$formula_X1, data = data)[, -1, drop = FALSE] 
    } else {
      X <- modelmatrixI(ff$formula_X1, data = data, XZ = FALSE)
    }
  }
  msd <- mean_sd_w(X, w = weights)
  if (novif) {
    return(list(msd = msd, ff_get = ff))
  }
  avg <- msd$mean
  sds <- msd$sd;
  sWn <- sqrt(weights)
  Zw <- model.matrix(ff$formula_Z, data) * sWn
  qrZ <- qr(Zw)
  if (XZ) {
    Xw <- qr.resid(qr(matrix(sWn)), X * sWn) 
  } else { 
    Xw <- qr.resid(qrZ, X * sWn)
  }
  qrX <- qr(Xw)
  diagXtX_inv <- diag(chol2inv(qrX$qr, size = qrX$rank))
  diagXtX_inv <- c(diagXtX_inv, rep(NA, length(sds)- length(diagXtX_inv)))
  VIF <- diagXtX_inv * sds ^ 2
  attr(qrX$qr, "scaled:center") <- c(msd$mean)
  attr(qrX$qr, "sd") <- c(msd$sd)
  meansdvif <- as.data.frame(t(rbind(avg, sds, VIF)))
  names(meansdvif) <- c("Avg", "SDS", "VIF")
  rr <- list(meansdvif = meansdvif, qrZ = qrZ, qrX = qrX, Zw = Zw, Xw = Xw)
  if (XZ) {
    names(rr)[c(3, 5)] <- c("qrXZ", "XZw")
  }
  return(rr)
}

# Hill number of order 2: N2
#' @noRd
#' @keywords internal
fN2 <- function(x) {
  x <- x / sum(x)
  1 / sum(x * x)
}
#

# from vegan 2.6-4 --------------------------------------------------------

#' @noRd
#' @keywords internal
centroids.cca <-  function(x, 
                           mf,
                           wt) {
  if (is.null(mf) || is.null(x)) {
    return(NULL)
  }
  facts <- sapply(mf, is.factor) | sapply(mf, is.character)
  if (!any(facts)) {
    return(NULL)
  }
  mf <- mf[, facts, drop = FALSE]
  ## Explicitly exclude NA as a level
  mf <- droplevels(mf, exclude = NA)
  if (missing(wt)) {
    wt <- rep(1, nrow(mf))
  }
  ind <- seq_len(nrow(mf))
  workhorse <- function(x, wt) {
    colSums(x * wt) / sum(wt)
  }
  ## As NA not a level, centroids only for non-NA levels of each factor
  tmp <- lapply(mf, function(fct) {
    tapply(ind, fct, function(i) workhorse(x[i, , drop = FALSE], wt[i]))
  })
  tmp <- lapply(tmp, function(z) sapply(z, rbind))
  pnam <- labels(tmp)
  out <- NULL
  if (ncol(x) == 1) {
    nm <- unlist(sapply(pnam,
                        function(nm) paste0(nm, names(tmp[[nm]]))),
                 use.names = FALSE)
    out <- matrix(unlist(tmp), nrow = 1, dimnames = list(NULL, nm))
  } else {
    for (i in seq_along(tmp)) {
      colnames(tmp[[i]]) <- paste0(pnam[i], colnames(tmp[[i]]))
      out <- cbind(out, tmp[[i]])
    }
  }
  out <- t(out)
  colnames(out) <- colnames(x)
  return(out)
}

#' @noRd
#' @keywords internal
modelmatrixI <- function(formula, 
                         data, 
                         XZ = TRUE){
  # model matrix with full identity contracts = full indicator coding
  ccontrasts <- function(x, data) {
    contrasts(data[[x]], contrasts = FALSE)
  }
  ff <- get_Z_X_XZ_formula(formula, data)
  if (XZ) {
    f <- ff$formula_XZ 
  } else {
    f <- ff$formula_X1
  }
  fcts <- c(ff$Condi_factor, ff$focal_factor)
  if (length(fcts)) {
    cntI <-  lapply(as.list(fcts), ccontrasts, data = data)
    names(cntI) <- fcts
    X <- model.matrix(f, contrasts.arg = cntI, data = data)
  } else {
    X <- model.matrix(ff$formula_XZ, data = data)
  }
  return(X[, -1, drop = FALSE])
}