#' @title Species-level Permutation Test in Double Constrained Correspondence Analysis (dc-CA)
#'
#' @description
#' \code{anova_species} performs a species-level permutation test of dc-CA.
#' The test uses residual predictor permutation (ter Braak 2022), which is robust
#' against differences in species total abundance in the \code{response} in \code{\link{dc_CA_vegan}} (ter Braak & te Beest, 2022)
#' The arguments of the function are similar to those of \code{\link[vegan]{anova.cca}},
#' but more restrictive.
#'
#' @param  object  an object from \code{\link{dc_CA_vegan}}.
#' @param permutations a list of control values for the permutations as
#'  returned by the function \code{\link[permute]{how}}, or
#'  the number of permutations required (default 999) or
#'  a permutation matrix where each row gives the permuted indices.
#'
#' @param  by character \code{"axis"} which sets the test statistic to the
#'  first eigenvalue of the dc-CA model.
#'  Default: \code{NULL} which set the test statistic to the inertia named \code{constraintsTE}
#'  in the inertia element of \code{\link{dc_CA_vegan}}). This is
#'  the environmentally constrained inertia explained by the traits (without trait covariates).
#'  (which is equal to the trait-constrained inertia
#'  explained by the environmental predictors (without covariates).)
#'  The default is quicker computationally as it avoid computation of an svd of permuted data sets.
#' @details
#' In \code{\link{anova_species}}, the first step extracts the species-niche centroids (SNC)
#' with respect to all dc-CA ordination axes
#' from an existing dc-CA analysis. These SNCs (more precisely, a rotated version thereof)
#' could have been obtained from a \code{\link[vegan]{cca}}
#' of the abundance data with the environmental variables as predictors, but this is computationally slower,
#' of course, than just extracting the unconstrained species scores from an existing dc-CA. The second step,
#' applies a weighted redundancy analysis of these SNCs with the traits as predictors. The second step is thus
#' a species-level analysis, but the final results (eigenvalues/ordination axes) are identical to those of the
#' analyses steps in \code{\link{dc_CA_vegan}}.The second step uses
#' published R-code for weighted redundancy analysis, which includes statistical significance
#' tests using residual predictor permutation (ter Braak, 2022).
#'
#' @return
#'  A list with two elements with names \code{table} and \code{eig}.
#'  The \code{table} is as from \code{\link[vegan]{anova.cca}} and \code{eig} gives the dc-CA eigenvalues,
#'  This output can be used for scripting forward selection of traits,
#'  similar to the forward selection of environmental
#'  variables in the demo \code{dune_select.r}.
#' @references
#' ter Braak, C.J.F. & te Beest, D.E. 2022. Testing environmental effects
#' on taxonomic composition with canonical correspondence analysis:
#' alternative permutation tests are not equal.
#' Environmental and Ecological Statistics. 29 (4), 849-868.
#' https://doi.org/10.1007/s10651-022-00545-4
#'
#' ter Braak, C.J.F. (2022) Predictor versus response permutation
#' for significance testing in weighted regression and redundancy analysis.
#' Journal of statistical computation and simulation, 92, 2041-2059.
#'  https://doi.org/10.1080/00949655.2021.2019256
#' @example demo/dune_test.R
#' @export



anova_species <- function(object, permutations = 999, by = NULL ){
# object dcca object
#
  if (is.null(by)) by <- "omnibus"
  if (is.na(pmatch(by, c("axis","omnibus")) )) stop(" set argument 'by' to 'axis' or 'NULL'")

  N <- nrow(object$data$dataTraits) #
  if (is.numeric(permutations)|| "how" %in% class(permutations) || is.matrix(permutations) ){
    if (is.numeric(permutations) && !is.matrix(permutations) ) permutations <- permute::how(nperm=permutations[1])
    else if (is.matrix(permutations) && !ncol(permutations)== N)
      stop(paste("Error:: each row of permutations should have ", N, " elements", sep =""))

  } else stop("argument permutations should be integer, matrix or specified  by permute::how(). ")




formulaEnv <- change_reponse(object$formulaEnv, "object$data$Y")
environment(formulaEnv)<- environment()
step1_sp <-vegan::cca(formulaEnv, data = object$data$dataEnv)

SNC_orthonormal_env <- vegan::scores(step1_sp, display= "species",
                                         scaling = "species", choices = 1:Rank_mod(step1_sp))
if (rownames(SNC_orthonormal_env)[1]=="col1") rownames(SNC_orthonormal_env) <- paste("Species", seq_len((nrow(object$data$dataTraits))),sep="")
#colSums(SNC_orthonormal_env^2 * object$weights$columns)/step1_sp$CCA$eig

# step 2 Perform a weighted RDAR(M^*~E): an RDA of M^* on the environmental variables-------------
#        using row weights R.

sWn <- sqrt(object$weights$columns)
Yw <-  SNC_orthonormal_env*sWn
Xw <- qr.X(get_QR(object$CCAonTraits,model = "CCA"))
qrZ <- get_QR(object$CCAonTraits, model = "pCCA")
if(is.null(qrZ)) Zw <- matrix(sWn) else  Zw<-  cbind(sWn, SVD(qr.X(qrZ)))


  # residual predictor permutation
  out_tes <- list()
  out_tes[[1]]  <- randperm_eX0sqrtw(Yw,Xw, Zw, sWn = sWn,  permutations= permutations, by = by, return= "all")

  if( by == "axis") {

  while (out_tes[[1]]$rank > length(out_tes) ) {
    Zw <- cbind(Zw,out_tes[[length(out_tes)]]$EigVector1)
    out_tes[[length(out_tes)+1]] <- randperm_eX0sqrtw(Yw,Xw, Zw,
                                  sWn = sWn, permutations= permutations, by = by, return = "all")
  }
  }
# what the env. variables explain of the trait-structured variation
  ss <-c(sapply(out_tes, function(x)x$ss[1]),out_tes[[length(out_tes)]]$ss[2])
  df <- #c(sapply(out_tes, function(x)x$ss[1]),out_tes[[length(out_tes)]]$ss[2])
    c(rep(1, length(ss)-1), out_tes[[length(out_tes)]]$df[2] )

  names(df) <- c(paste("dcCA", seq_along(out_tes), sep = ""),"Residual")
  fraqExplained <-c(sapply(out_tes, function(x)x$ss[1])/sum(out_tes[[1]]$ss),NA)

  F0 <-c(sapply(out_tes, function(x)x$F0[1]),NA)
  F.perm <- out_tes[[1]]$Fval
  if (length(out_tes)>1){
    for (k in seq_along(out_tes)[-1]) F.perm <- cbind(F.perm, out_tes[[k]]$Fval )
  }


  p_val_axes1 <- c(cummax(sapply(out_tes, function(x)x$pval[1])),NA)
  eig <- out_tes[[1]]$eig
  axsig_dcCA_species <- data.frame(df =df, `ChiSquare`= ss,R2 = fraqExplained, `F`  = F0,`Pr(>F)` = p_val_axes1)
  names(axsig_dcCA_species)[5]<- "Pr(>F)"


  #data.frame(f_species0$df, f_species0$ss, c(f_species0$F0, NA), c(f_species0$pval, NA))
  #colnames(table) <- c("Df", varname, "F", "Pr(>F)")
  object1 <- paste("Model:", c(object$call),"\n")
  head <- paste0("Species-level permutation test using dc-CA\n",
                 object1,
                 "Residualized predictor permutation\n",
                 howHead(attr(out_tes[[1]],"control") ))

  f_species <-structure(axsig_dcCA_species, heading = head, #Random.seed = seed,
                        control = attr(out_tes[[1]],"control"),
                        Random.seed =   attr(out_tes[[1]],"seed"),
                        control = attr(out_tes[[1]],"control"),
                        F.perm = F.perm,
                        class = c("anova.cca", "anova", "data.frame"))

  result <- list(table = f_species, eig = eig)
return(result)
}


# functions from JSCS 2022------------------------------------------------------

# function dummysvd is to avoid calculation of the eigen vector/
#                                                  eigen value test statistic
#           set: mysvd <- dummysvd # with eigenvalue test statistic
# otherwise set  mysvd <- svd      # witout eigenvalue test statistic
dummysvd <- function(Y){list(d = 1,
                             U = matrix(1,nrow=nrow(Y), ncol=1),
                             V = matrix(1,nrow=ncol(Y), ncol=1) )}


#for permutation.type = X = X1
randperm_eX0sqrtw <- function(Y,X, Z = matrix(1, nrow = nrow(Y),ncol =1), by= NULL,
                              sWn = rep(1,nrow(Y)), permutations= permute::how(nperm=999),
                              return = "pval"){

  # for getting ordinary X-permutation, i.e. ordinary Collins-Dekker (residualized predictor permutation),
  #    permuting the X-residuals of the original weighted analysis,
  #    in the calling function
  # Y X Z should have been transformed to LS in the calling function by
  #  sWn = sqrt(w/sum(w)) of the calling function
  #  so that this function does not need weights.
  #  except for calculating the original residual X given Z
  #
  EPS <- sqrt(.Machine$double.eps) # for permutation P-values
  # preparations and data value
  if (by =="axis") mysvd<-svd else mysvd <- dummysvd
  N <- nrow(Y)
  if (is.matrix(permutations)){
    # matrix: check that it *strictly* integer
    if(!is.integer(permutations) && !all(permutations == round(permutations)))
      stop("permutation matrix must be strictly integers: use round()")
    perm.mat <- permutations
  } else if ("how" %in% class(permutations)){
    #perm.mat creation
    perm.mat <- permute::shuffleSet(N, control = permutations)
  }
  nrepet <- nrow(perm.mat)
  Z_orth <- SVD(Z)
  # Y-residuals from Y~Z under null model
  eY <- unweighted_lm_Orthnorm(Y, Z_orth)
  # orthogonalize X with respect to Z giving eX
  # eX = X_orth = residuals of X from X~Z
  eX <- unweighted_lm_Orthnorm(X, Z_orth)
  eXw <- eX/sWn # the X-residuals of the original weighted analysis
  # step 2: ss(X)
  Yfit_X <- unweighted_lm_pq_fit(eY,eX)
  ssX <- sum(Yfit_X^2)
  svd_Yfit_X <- svd(Yfit_X)
  ssX_eig1 <- svd_Yfit_X$d[1]^2
  EigVector1 <- svd_Yfit_X$u[, 1, drop= FALSE]
  rank <- sum(svd_Yfit_X$d/sum(svd_Yfit_X$d) > 1.e-3)
  sstot <- sum(eY^2)
  p_eX <- qr(eX)$rank
  df_cor <- (nrow(eY)- ncol(Z_orth)- p_eX)/p_eX  #(N-nz-nx-1)/nx
  F0 <- ssX / (sstot - ssX) * df_cor
  F0_eig1 <- ssX_eig1 / (sstot -  ssX_eig1) * df_cor
  # end preparations and data value
  ssX_perm <- numeric(nrepet)
  ssX_eig1_perm <- numeric(nrepet)
  #i=1
  #
  for (i in seq_len(nrepet)){
    # permute residual eX
    i_perm  <- perm.mat[i,]
    eX_perm <- eXw[i_perm, , drop = FALSE]*sWn
    # orthogonalize permuted orthogonaled-X
    #  (i.e. eX which is permuted to give eX_perm) with respect to Z
    #  giving eX_perm_orth_to_Z
    # residuals of X_orth_perm on ~Z :
    eX_perm_orth_to_Z <- unweighted_lm_Orthnorm(eX_perm, Z_orth)
    # regress the Y-residuals w.r.t. Z (i.e. eY) on eX_perm_orth_to_Z
    # to determine the sum of squares due to eX_perm in the regression lm(Y~ Z + eX_perm), i.e.
    # for one response variable:
    # my_lm <- lm(Y~ Z + eX_perm, weights = Wn); anova(my_lm)['eX_perm',2]
    Yfit_permX <- unweighted_lm_pq_fit(eY,eX_perm_orth_to_Z) #
    ssX_perm[i] <- sum(Yfit_permX^2)
    ssX_eig1_perm[i] <- mysvd(Yfit_permX)$d[1]^2
  }
  if (by =="axis")  {
    ssX_perm<- ssX_eig1_perm
    F0 <- F0_eig1
    }

  rss_perm <- sstot-ssX_perm
  Fval <- ssX_perm /rss_perm * df_cor
  isna.r <-  sum(is.na(Fval))
  pval <- (sum(Fval >= (F0 - EPS), na.rm = TRUE) + 1)  / (nrepet- isna.r  + 1)
  attr(pval, "test") <- by
  if (return == "pval") res <- pval else {

    eig <- svd_Yfit_X$d
    eig <- (eig[eig>EPS])^2
    res <- list(pval = pval, Fval = Fval, F0 = F0,
              ss = c(Model= ssX, Residual = sstot-ssX),
              df = c(Model = qr(eX)$rank, Residual = nrow(eY)- ncol(Z_orth)- p_eX ),
              rank = rank,
              eig = eig, EigVector1 = EigVector1)
    if (is.matrix(permutations)){
            attr(perm.mat, "control") <-
              structure(list(within=list(type="supplied matrix"),
                             nperm = nrow(perm.mat)), class = "how")
    }
    attr(res, "seed") <- attr(perm.mat,"seed")
    attr(res, "control") <-  attr(perm.mat, "control")
  }

return(res)
}

# unweighted least-squares (OLS) functions -------------------------------------------------

unweighted_lm_pq_fit <- function(Y, X){
  # multivariate multiple regression of Y on X
  # using qr decomposition
  # value Y_fit
  Y_fit <- qr.fitted( qr(X),Y)
  return(Y_fit)
}

unweighted_lm_Orthnorm_fit <- function(Y, X=numeric(0)){
  # multivariate multiple regression of Y on orthonormal X
  # value Y_residual
  beta <-  t(X) %*% Y
  Yfit <- X %*% beta
  return(Yfit)
}

unweighted_lm_Orthnorm <- function(Y, X=numeric(0)){
  # multivariate multiple regression of Y on orthonormal X
  # value Y_residual
  beta <-  t(X) %*% Y
  Y <- Y - X %*% beta
  return(Y)
}



SVD <- function(Y){
  svdY <- svd(Y)
  Ustar <- svdY$u
  id <- which(svdY$d > 1.e-6)
  return(Ustar[,id, drop = FALSE])
}

# end OLS versions -------------------------------------------------


# code adapted from vegan 2.6-4 ---------------------------------------------------

### Make a compact summary of permutations. This copies Gav Simpson's
### permute:::print.how, but only displays non-default choices in how().
`howHead` <- function(x, ...)
{
  ## print nothing is this not 'how'
  if (is.null(x) || !inherits(x, "how"))
    return()
  ## collect header
  head <- NULL
  ## blocks
  if (!is.null(permute::getBlocks(x)))
    head <- paste0(head, paste("Blocks: ", x$blocks.name, "\n"))
  ## plots
  plotStr <- permute::getStrata(x, which = "plots")
  if (!is.null(plotStr)) {
    plots <- permute::getPlots(x)
    ptype <- permute::getType(x, which = "plots")
    head <- paste0(head, paste0("Plots: ", plots$plots.name, ", "))
    head <- paste0(head, paste("plot permutation:", ptype))
    if(permute::getMirror(x, which = "plots"))
      head <- paste(head, "mirrored")
    if (isTRUE(all.equal(ptype, "grid"))) {
      nr <- permute::getRow(x, which = "plots")
      nc <- permute::getCol(x, which = "plots")
      head <- paste0(head, sprintf(ngettext(nr, " %d row", " %d rows"),
                                   nr))
      head <- paste0(head, sprintf(ngettext(nc, " %d column",
                                            " %d columns"), nc))
    }
    head <- paste0(head, "\n")
  }
  ## the fine level (within plots if any)
  type <- permute::getType(x, which = "within")
  head <- paste0(head, "Permutation: ", type)
  if (isTRUE(type %in% c("series", "grid"))) {
    if(permute::getMirror(x, which = "within"))
      head <- paste(head, "mirrored")
    if(permute::getConstant(x))
      head <- paste0(head, " constant permutation within each Plot")
  }
  if (isTRUE(all.equal(type, "grid"))) {
    nr <- permute::getRow(x, which = "within")
    nc <- permute::getCol(x, which = "within")
    head <- paste0(head, sprintf(ngettext(nr, " %d row", " %d rows"),
                                 nr))
    head <- paste0(head, sprintf(ngettext(nc, " %d column",
                                          " %d columns"), nc))
  }
  paste0(head, "\nNumber of permutations: ", permute::getNperm(x),  "\n")
}

