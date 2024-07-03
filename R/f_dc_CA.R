#' @noRd
#' @keywords internal
change_reponse  <- function(f, 
                            response, 
                            data = NULL){ # used in dc_CA and dc_CA
  # response : character
  ft <- as.character(f)
  if (any(ft %in% ".")) {
    # change the dot using  data
    if (is.null(data)) {
      stop("data should be present if formula contains a period(.).")
    }
    if (ft[2] == ".") { 
      ft <- c("~", paste(names(data), collapse = "+")) 
    } else {
      stop("specify the formulas with explicit remaining predictors")
    }
  }
  if (!is.character(response)) { 
    stop("response must be character")
  } else if (length(response) > 1) { 
    stop("response must be of length one")
  }
  id <- which(ft %in% "~")
  if (!length(id)) {
    stop("specify a formula with a '~' ")
  } else if (length(ft) > 3) {
    warning("as.character(formula) has more than three elements, ", 
            "only its last element is included")
  }
  ft2 <- paste0(response," ", ft[id], ft[length(ft)], collapse = " ")
  f2 <- as.formula(ft2)
  return(f2)
}

#' @noRd
#' @keywords internal
calculate_b_se_tval <- function(X_or_qr_decomp_of_X, 
                                y, 
                                w = NULL, 
                                scale2 = 0, 
                                name = "SNC", 
                                fitted_only = FALSE) {
  # specify  X_or_qr_decomp_of_X is an (yet unweighted to-be weighted) matrix or
  # the qr_decomp of the weighted X matrix
  # y is the (yet unweighted to-be weighted) response vector or matrix
  # wrda TRUE arguments come from function wrda, If from vegan FALSE
  if (is.null(w)) {
    if (is.matrix(X_or_qr_decomp_of_X)) {
      w <- rep(1, nrow(X_or_qr_decomp_of_X))
    } else {
      w <- rep(1, nrow(X_or_qr_decomp_of_X$qr))
    }
  }
  w <- w / sum(w) 
  sqrtw <- sqrt(w)
  y_weightedw <- as.matrix(y) * sqrtw
  TSS <- colSums(y_weightedw ^ 2)
  # Perform QR decomposition on the weighted design matrix
  if (!is.qr(X_or_qr_decomp_of_X)) {
    X_weighted <- X_or_qr_decomp_of_X * sqrtw
    QR <- qr(X_weighted)
  } else {
    QR <- X_or_qr_decomp_of_X
  }
  # Compute estimated regression coefficients
  beta_hat <- qr.coef(QR, y_weightedw)
  # Calculate fitted values
  fitted_valuesw <- qr.fitted(QR, y_weightedw)
  # Calculate residuals
  residuals <- y_weightedw - fitted_valuesw
  fitted_values <- fitted_valuesw / sqrtw
  if (fitted_only) {
    out1 <- fitted_values
  } else {
    # Compute residual sum of squares (RSS)
    RSS <- colSums(residuals ^ 2)
    # Estimate variance of the errors
    n <- length(w)
    p <- QR$rank
    sigma_hat_sq <- RSS / (n - p - 1 ) 
    # Calculate variance-covariance matrix of the estimated regression coefficients
    diagXtX_inv <- diag(chol2inv(QR$qr, size = QR$rank))
    se <- matrix(nrow = nrow(beta_hat), ncol = length(sigma_hat_sq))
    for (i in seq_along(sigma_hat_sq)) {
      var_covar_matrix <- diagXtX_inv * sigma_hat_sq[i]
      # Calculate standard errors
      se[, i] <- sqrt(var_covar_matrix)
    }
    TSSfit <- colSums(fitted_valuesw ^ 2)
    if (scale2) {
      sqrtTSS <- sqrt(TSSfit / scale2)
      if (length(TSS) > 1) {
        sTSS <- diag(1 / sqrtTSS) 
      } else { 
        sTSS <- diag(1) * (1 / sqrtTSS)
      }
      fitted_values <- fitted_values %*% sTSS
      y <- y %*% sTSS
      beta_hat <- beta_hat %*% sTSS
      se <- se %*% sTSS
    }
    if (length(TSS) > 1) {
      colnames(y) <- paste0(name, seq_len(ncol(y)))
      colnames(fitted_values) <- paste0(name, "_lc", seq_len(ncol(y)))
    }
    avg <- attr(QR$qr, which = "scaled:center")
    sds <- sqrt(colSums(qr.X(QR) ^ 2))
    tval <- beta_hat / se
    VIF <- diagXtX_inv * sds ^ 2
    beta_stan <- beta_hat * sds
    colnames(beta_stan) <- paste0("Regr", seq_len(ncol(beta_stan)))
    colnames(tval) <- paste0("tval", seq_len(ncol(tval)))
    b_se <- data.frame(beta = beta_hat, se = se)
    coef_normed <- cbind(Avg = avg, SDS = sds, VIF = VIF, beta_stan, tval)
    attr(coef_normed, "meaning") <- 
      "mean, sd, VIF, standardized regression coefficients and their t-ratio"
    out1 <- list(fitted = fitted_values, y = y, coef_normed = coef_normed, 
                 b_se = b_se,  R2 = 1 - RSS / TSS)
  }
  return(out1)
}
