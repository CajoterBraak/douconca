f_inertia <- function(object){
  # function f_inertia uses vegan 2.6-4 internal structure
  # object: a dccav object, results of dc_CA_vegan
  # value:  a matrix (currently with 1 column) with in first column the inertias
  inertia <- cbind(c(total=           object$CCAonTraits$tot.chi,
                     conditionT=      object$CCAonTraits$pCCA$tot.chi,
                     traits_explain=  object$CCAonTraits$CCA$tot.chi,
                     conditionE =     object$RDAonEnv$pCCA$tot.chi,
                     constraintsTE =   object$RDAonEnv$CCA$tot.chi ))
  colnames(inertia)<- "weigthed variance"
  expla <- c("total inertia","inertia of the trait condition", "trait-constrained inertia",
             "trait-constrained inertia explained by the condition in formulaEnv",
             "trait-constrained inertia explained by the predictors in formulaEnv")
  names(expla) <- c("total","conditionT","traits_explain","conditionE","constraintsTE")
  attr(inertia, which = "meaning") <-  matrix( expla[rownames(inertia)], ncol=1,
                                               dimnames = list(rownames(inertia),"meaning"))
  return(inertia)
}

get_QR <- function(object, model = "CCA"){
  # function get_QR uses vegan 2.6-4 internal structure
  # gets the qr decompostion of object
  # model = "CCA" or "pCCA"
  if (model=="CCA")  QR<- object$CCA$QR else if (model == "pCCA") QR <- object$pCCA$QR  else stop("model not supported")
  return(QR)
}




