#' @title Community-level Permutation Test in Double Constrained Correspondence Analysis (dc-CA)
#'
#' @description
#' \code{anova_sites} performs a Community-level permutation test of dc-CA.
#' The test uses residual predictor permutation (ter Braak 2022), which is robust
#' against differences in sites total abundance in the \code{response} in \code{\link{dc_CA_vegan}} (ter Braak & te Beest, 2022)
#' The arguments of the function are similar to those of \code{\link[vegan]{anova.cca}}, but more restricted.
#' With equal site-totals as in \code{\link{dc_CA_vegan}}, \code{anova(object$RDAonEnv)} is much faster.
#
#' @param  object  an object from \code{\link{dc_CA_vegan}}.
#' @param permutations a list of control values for the permutations as
#'  returned by the function \code{\link[permute]{how}}, or
#'  the number of permutations required (default 999), or
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
#' The algorithm is two-step. The first step is a \code{\link[vegan]{cca}}
#' of the transposed \code{response} on to the traits using
#' \code{formulaTraits}. The second is a weighted redundancy analysis of the community weighted means (CWM) of
#' orthonormalized traits againt the environemtal variables using \code{formulaEnv} and
#' published R-code for weighted redundancy analysis, which includes statistical significance
#' tests using residual predictor permutation (ter Braak, 2022).
#'
#' @return
#'  A list with two elements with names \code{table} and \code{eig}.
#'  The \code{table} is as from \code{\link[vegan]{anova.cca}} and \code{eig} gives the dc-CA eigenvalues,
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



anova_sites <- function(object, permutations =999, by = NULL){
# object dcca object;
#
  if (is.null(by)) by <- "omnibus"
  if (is.na(pmatch(by, c("axis","omnibus")) )) stop(" set argument 'by' to 'axis' or 'NULL'")

  N <- nrow(object$data$dataEnv) #
  if (is.numeric(permutations)|| "how" %in% class(permutations) || is.matrix(permutations) ){
      if (is.numeric(permutations) && !is.matrix(permutations) ) permutations <- permute::how(nperm=permutations[1])
      else if (is.matrix(permutations) && !ncol(permutations)== N)
        stop(paste("Error:: each row of permutations should have ", N, " elements", sep =""))

  } else stop("argument permutations should be integer, matrix or specified  by permute::how(). ")


# start of new dc-ca ------------------------------------------------------

  step1 <- object$CCAonTraits
  out1 <- object[c("CCAonTraits", "formulaTraits","data", "weights")]

  n <- nrow(out1$data$Y)
  CWMs_orthonormal_traits <- vegan::scores(step1, display= "species",
                                           scaling = "species", choices = 1:Rank_mod(step1)) #* sqrt((n-1)/n)
  if (rownames(CWMs_orthonormal_traits)[1]=="col1") rownames(CWMs_orthonormal_traits) <- paste("Sam", seq_len((nrow(out1$data$dataEnv))),sep="")

# step 2 Perform a weighted RDAR(M^*~E): an RDA of M^* on the environmental variables-------------
#        using row weights R.

sWn <- sqrt(object$weights$rows)
Yw <-  CWMs_orthonormal_traits*sWn
Xw <- qr.X(get_QR(object$RDAonEnv,model = "CCA"))
qrZ <- get_QR(object$RDAonEnv, model = "pCCA")
if(is.null(qrZ)) Zw <- matrix(sWn) else  Zw<-  cbind(sWn, SVD(qr.X(qrZ)))


  # residual predictor permutation
  out_tes <- list()
  out_tes[[1]]  <- randperm_eX0sqrtw(Yw,Xw, Zw, sWn = sWn, permutations= permutations, by = by, return = "all")

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
  axsig_dcCA_sites <- data.frame(df =df, `ChiSquare`= ss,R2 = fraqExplained, `F`  = F0,`Pr(>F)` = p_val_axes1)
  names(axsig_dcCA_sites)[5]<- "Pr(>F)"

  object1 <- paste("Model:", c(object$call),"\n")

  head <- paste0("Community-level permutation test using dc-CA\n",
                 object1,
                 "Residualized predictor permutation\n",
                 howHead(attr(out_tes[[1]],"control") ))

  f_sites <-structure(axsig_dcCA_sites, heading = head, #Random.seed = seed,
                        control = attr(out_tes[[1]],"control"),
                        Random.seed =   attr(out_tes[[1]],"seed"),
                        control = attr(out_tes[[1]],"control"),
                        F.perm = F.perm,
                        class = c("anova.cca", "anova", "data.frame"))

  result <- list(table = f_sites, eig = eig)

return(result)
}
