# utility functions for f_trait_axes and f_env_axes in print.dcca and scores.dcca--------------


wcor <- function(X, Y=X, w = rep(1,nrow(X))){
  # weighted correlation between matrix X and Y
  w <- w/sum(w)
  Xstd <- standardize_w(X, w)
  Ystd <- standardize_w(Y, w)
  t(Xstd) %*% diag(w) %*% Ystd
}

mean_w <- function(X,w = rep(1/nrow(X),nrow(X))){t(w/sum(w))%*% X}
center_w <- function(X,w = rep(1/nrow(X),nrow(X))){ X - rep(1,length(w))%*%t(w)%*% X }

#standardize_w <- function(X,w = rep(1/nrow(X),nrow(X)), wsvd = FALSE){
standardize_w <- function(X,w = rep(1/nrow(X),nrow(X))){
  # NB requires w to be have sum 1
  ones <- rep(1,length(w))
  Xc <- X - ones %*% t(w)%*% X
  Xstd <- Xc / ones%*%sqrt(t(ones)%*%(Xc*Xc*w))
  #if (wsvd) Xstd <- wSVD(Xstd, w)
  return(Xstd)
}

mean_sd_w <- function(X,w = rep(1/nrow(X),nrow(X))){
  # NB requires w to be have sum 1
  ones <- rep(1,length(w))
  mean_w <- t(w/sum(w))%*% X
  Xc <- X - ones %*% mean_w
  sd_w <-  sqrt(t(ones)%*%(Xc*Xc*w))
  #if (wsvd) Xstd <- wSVD(Xstd, w)
  return(list(mean = mean_w, sd = sd_w))
}

msdvif <- function(formula = NULL, data, weights, XZ = FALSE){
  # calc mean variance and vif from for X given Z or XZ with qr of X|Z or of centered XZ
  if (is.null(formula)) {f <- ~. } else {f <- formula}
  ff <- get_Z_X_XZ_formula(f,data)
  if (XZ)X <- model.matrix(ff$formula_XZ, data = data)[,-1, drop = FALSE] else
         X <- model.matrix(ff$formula_X1, data= data)[,-1, drop = FALSE]
  msd <- mean_sd_w(X,w= weights)
  avg = msd$mean ; sds = msd$sd;
  sWn <- sqrt(weights)
  Zw <- model.matrix(ff$formula_Z, data)*sWn
  qrZ <- qr(Zw)
  if (XZ) Xw <- qr.resid(qr(matrix(sWn)), X*sWn) else Xw <- qr.resid(qrZ, X*sWn)
  qrX <- qr(Xw)
  diagXtX_inv <-   diag(chol2inv(qrX$qr, size = qrX$rank))
  VIF <- diagXtX_inv*sds^2
  attr(qrX$qr, "scaled:center") <- c(msd$mean)
  attr(qrX$qr, "sd") <- c(msd$sd)
  meansdvif <- as.data.frame(t(rbind(avg,  sds,  VIF)))
  names(meansdvif) <- c("Avg", "SDS","VIF")
  rr <- list(meansdvif=meansdvif, qrZ= qrZ, qrX = qrX, Zw= Zw, Xw =Xw)
  if (XZ)names(rr)[c(3,5)]<- c("qrXZ", "XZw")
  return(rr)
}

# Hill number of order 2: N2
fN2 <- function(x){x <- x/sum(x); 1/sum(x*x)}
#

# from vegan 2.6-4 --------------------------------------------------------

centroids.cca <-  function(x, mf, wt)
  {
    if (is.null(mf) || is.null(x))
      return(NULL)
    facts <- sapply(mf, is.factor) | sapply(mf, is.character)
    if (!any(facts))
      return(NULL)
    mf <- mf[, facts, drop = FALSE]
    ## Explicitly exclude NA as a level
    mf <- droplevels(mf, exclude = NA)
    if (missing(wt))
      wt <- rep(1, nrow(mf))
    ind <- seq_len(nrow(mf))
    workhorse <- function(x, wt)
      colSums(x * wt) / sum(wt)
    ## As NA not a level, centroids only for non-NA levels of each factor
    tmp <- lapply(mf, function(fct)
      tapply(ind, fct, function(i) workhorse(x[i,, drop=FALSE], wt[i])))
    tmp <- lapply(tmp, function(z) sapply(z, rbind))
    pnam <- labels(tmp)
    out <- NULL
    if (ncol(x) == 1) {
      nm <- unlist(sapply(pnam,
                          function(nm) paste(nm, names(tmp[[nm]]), sep="")),
                   use.names=FALSE)
      out <- matrix(unlist(tmp), nrow=1, dimnames = list(NULL, nm))
    } else {
      for (i in seq_along(tmp)) {
        colnames(tmp[[i]]) <- paste(pnam[i], colnames(tmp[[i]]),
                                    sep = "")
        out <- cbind(out, tmp[[i]])
      }
    }
    out <- t(out)
    colnames(out) <- colnames(x)
    out
  }

