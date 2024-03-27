#' @title Get constraining and conditional variables from an rda or cca ordination
#'
#' @description
#' \code{get_focal_and_conditioning_factors} derives the focal constraining variables(s) and (if present) the
#' conditioning variables(s) from a \code{\link[vegan]{rda}}
#' or \code{\link[vegan]{cca}} ordination specified via a formula for use in a PRC
#' (with \code{\link{doPRC}} or \code{\link{PRC_scores}} or in the dcCA library for dc-CA scores with scores.dccav).
#'
#' @param  object  a result of \code{\link[vegan]{rda}}, \code{\link[vegan]{cca}} specified via a formula (S3 method for class 'formula'),
#' or a result of \code{\link{doPRC}} or \code{\link{PRC_scores}}.
#' @param factors_only logical to include factors only. If the formula contains one or more interactions,
#' only the factors appearing in the interaction are returned if not appearing in the \code{Condition}.
#' Default \code{TRUE} which makes sense for PRC.
#' @return  A list with element \code{`focal factor`} and \code{condition}
#' @example demo/PRC_pyrifos_bk.R
#' @references
#' ter Braak C.J.F. & Šmilauer P. (2018): Canoco reference manual and user's guide:
#'  software for ordination, version 5.1x. Microcomputer Power, Ithaca, USA, 536 pp.
#'  (http::www.canoco5.com)
#' @seealso \code{\link{doPRC}}, \code{\link{PRC_scores}}
#' @noRd
get_focal_and_conditioning_factors <- function(object, factors_only = TRUE){
  #get_focal_and_conditioning_factors  from an cca-object
  tl <- attr(object$terms, "term.labels")
  if (is.null(tl)){
    out <- list(NULL, NULL)
    names(out) <-  c("focal factor", "condition")
    warning("Specify rda or cca via a formula, so that focal factor and condition can be found.")
    return(out)
  }
  # find all variable names by deleting from tl Condition and all interactions
  # first delete Condition
  idC <- pmatch("Condition(" , tl)
  if (is.na(idC)) { Condi_nams <- NULL} else {
    # get names of the Condition
    condi <- tl[idC]
    condi <- attr(stats::terms(stats::as.formula(paste("~",substr(condi, 11, nchar(condi)-1), sep = ""))), "term.labels")
    tl2 <- strsplit(condi, split = ":", fixed = TRUE)
    id_interactions <- which (sapply(tl2, length) >1)
    if (length(id_interactions)) Condi_nams <- unlist(tl2[-id_interactions]) else Condi_nams <- unlist(tl2)
    tl <- tl[-idC]

  }
  # Condi_nams set, now work on reduced rank constraining part
  # As it is a PRC look for an interaction
  tl2 <- strsplit(tl, split = ":", fixed = TRUE)
  if (!factors_only){
    #focal_nams <- unlist(tl2[which (sapply(tl2, length) ==1)])
    focal_nams <- unique(c(unlist(tl2[which (sapply(tl2, length) ==1)]),
             unlist(tl2[which (sapply(tl2, length) >1)])))
    if (!is.na(idC))focal_nams<- focal_nams[!focal_nams %in% condi]
  } else {

    # which interaction has max length, or 1 of no interaction
    id_interactions <- which.max (sapply(tl2, length))
    if (id_interactions==1 && length(tl2[[id_interactions]])==1){
      # no interactions
      interaction <- NULL;
      #focal_nams <- setdiff(unlist(tl2), Condi_nams)
      id <- unlist(tl2) %in% Condi_nams
      focal_nams <- unlist(tl2)[!id]
      if (factors_only){
        focal_nams <- focal_nams[focal_nams %in% names(object$terminfo$xlev)]
      }
     # if (length(focal_nams)==1) {Condi_nams <- focal_nams} else Condi_nams <-NULL
    } else {
      interaction <-  unlist(tl2[id_interactions])
      id <- which(interaction %in% Condi_nams)
      if (length(id)) {
        focal_nams <- interaction[-id]; Condi_nams <- interaction[id]
       # if (factors_only){
        # there may be numeric variables in focal_nams (e.g. if condition is a factor but numeric in treatment)
        # delete those
        # factors are names(object$terminfo$xlev)
        focal_nams <- focal_nams[focal_nams %in% names(object$terminfo$xlev)]
       # }
      }else {
        # unusual case:
        focal_nams <- interaction
        warning("interaction found in constraints, but none with conditional term")
      }
    }
  }
  out <- list(focal_nams, Condi_nams)
  names(out) <-  c("focal factor", "condition")
  attr(out, "interaction") <- !is.null(interaction)
  return(out)
}

