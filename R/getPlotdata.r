#'  Utility function: extracting data from a \code{\link{dc_CA_vegan}} object for plotting
#'  a single axis.
#' @description
#'  \code{getPlotdata} extracts data from a \code{\link{dc_CA_vegan}} object for plotting the
#' CWMs and SNCs of a single axis.
#' @param object results from \code{\link{dc_CA_vegan}} of class \code{dccav}.
#' @param axis the axis number to get (default 1).
#' @param envfactor name of row factor to display as color and lines in the CWM plot (default \code{NULL}).
#' The default extracts the factor from the environmental model.
#' If set to \code{NA}, no additional coloring and lines are displayed.
#' @param traitfactor name of column factor to display as color and lines in the SNC plot (default \code{NULL}).
#' The default extracts the factor from the trait model.
#' If set to \code{NA}, no additional coloring and lines are displayed.
#' @param newnames a list  with two elements:
#' names for traits and for environmental variables, default \code{NULL} for
#' names derived from the result of \code{\link{scores.dccav}} with \code{tidy=TRUE}.
# @param size.centroids size of centroid labels and points
#' @param remove.centroids logical to remove any centroids from the plot data (default \code{FALSE}).
#' Can be a two-vector, \emph{e.g.} \code{c(TRUE, FALSE)} to remove only the environmental centroids.
#' @param facet logical. Default \code{TRUE} for CWMs and SNCs plots in separate panels.
#' If \code{FALSE},this parameter changes the position of the environmental centroid names
#' (from left to right).
#' @example demo/dune_plot_dcCA.R
#' @export
#'
#'
getPlotdata <- function(object, axis=1, envfactor=NULL, traitfactor=NULL,newnames = NULL,
                      facet= TRUE,       remove.centroids=FALSE){
  size.centroids = 1

  # envfactor <- NULL # name or NA
  # traitfactor <- NULL # name or  NA
  # remove.centroids <- c(FALSE,FALSE) # todo

  # getPlotdata function

  if (length(remove.centroids)==1) remove.centroids <- c(remove.centroids ,remove.centroids)

  traitINcondition <- envINcondition <- FALSE
  if (is.null(envfactor)){
    ff <- get_focal_and_conditioning_factors(object$RDAonEnv, factors_only = TRUE)#'fo'
    if (is.null(ff$condition)){envfactor <- ff$`focal factor`[1]; envINcondition <- FALSE }else{
      envfactor <- ff$condition[1]; envINcondition <- TRUE
    }
  }

  if (is.null(traitfactor)){
    ff <- get_focal_and_conditioning_factors(object$CCAonTraits, factors_only = TRUE)#'fo'
    if (is.null(ff$condition)){traitfactor <- ff$`focal factor`[1]; traitINcondition <- FALSE }else{
      traitfactor <- ff$condition[1]; traitINcondition <- TRUE
    }
  }

  # end of set env and traitfactor

  mod_scores <- scores(object,choices=axis, tidy= TRUE, scaling = "symmetric")
  newname.list <- setnames(mod_scores, newnames = newnames)

  idTFc <- mod_scores$score %in% c("constraints_sites", "constraints_species", "centroids", "centroids_traits")
  idTFu <- mod_scores$score %in% c("sites", "species", "centroids", "centroids_traits")

  con_scores <- mod_scores[idTFc,]
  uncon_scores <- mod_scores[idTFu,]

  uncon_scores$score <- factor(uncon_scores$score)
  colnames(uncon_scores)[1]<- "CWM-SNC"

  trait_env_scores <- mod_scores[!(idTFc|idTFu),]


  scorepair <- cbind(con_scores[,-4],uncon_scores[,c(1,4),drop= FALSE])
  scorepair$score <- factor(scorepair$score)
  names(scorepair)[1]<-"dcCA1"
  # this removes centroids if traitfactor = NA
  if (remove.centroids[1] ) scorepair <- scorepair[!scorepair$score%in%c("centroids_traits"),]
  if (remove.centroids[2] ) scorepair <-  scorepair[!scorepair$score%in%c("centroids"),]

  if (!is.na(envfactor)|| is.character(envfactor)     ) {
    if (is.na(envfactor)) envfactor1 <- envfactor else if(envfactor %in% names(object$data$dataEnv) ) envfactor1 <- object$data$dataEnv[[envfactor]]
    else {stop(paste(envfactor, " must be in ", paste0(names(object$data$dataEnv),collapse = ",")))}
    } else envfactor1 <- envfactor


  if (!is.na(traitfactor)|| is.character(traitfactor) ) {
    if (is.na(traitfactor)) traitfactor1 <- traitfactor else if(traitfactor %in% names(object$data$dataTraits) ) traitfactor1 <- object$data$dataTraits[[traitfactor]]
    else {stop(paste(traitfactor, " must be in ", paste0(names(object$data$dataTraits),collapse = ",")))}
    } else traitfactor1<- traitfactor


  if (facet) {
    scorepair$dcCA1[scorepair$score%in% c("centroids","centroids_traits")] <-  min(scorepair$dcCA1,na.rm = TRUE)
  } else {
    scorepair$dcCA1[scorepair$score%in% c("centroids_traits")] <-  min(scorepair$dcCA1,na.rm = TRUE)
    scorepair$dcCA1[scorepair$score%in% c("centroids")] <-  max(scorepair$dcCA1,na.rm = TRUE)
  }



  #### could be integrated in below....

  if (nlevels(scorepair$score)==4) {
    typenams <- c("site", "taxon","centroid","centroid")
    scorepair$score <- factor(scorepair$score, levels = levels(scorepair$score)[c(3,4,1,2)])
  }else
    if (nlevels(scorepair$score)>2){
      typenams <- c("site", "taxon","centroid")
      scorepair$score <- factor(scorepair$score, levels = levels(scorepair$score)[c(2,3,1)])
    }else typenams <- c("site", "taxon")


  if (nlevels(scorepair$score)==4) typenams3 <- rep(c("CWM of trait composite", "SNC of dc-CA axis"),2) else
    if (nlevels(scorepair$score)==2) typenams3 <- c("CWM of trait composite", "SNC of dc-CA axis") else
      if (nlevels(scorepair$score)==3){
        if ("centroids_traits" %in% levels(scorepair$score)) typenams3 <- c("CWM of trait composite", rep("SNC of dc-CA axis",2))
        else typenams3 <- c("CWM of trait composite", "SNC of dc-CA axis","CWM of trait composite")
      }

  scorepair$type <-  typenams3[scorepair$score]
  #unique(scorepair$type)# "CWM of trait composite" "SNC of dc-CA axis"

  #scorepair$centroid <- factor(levels(scorepair$score)[scorepair$score])
  #scorepair$centroid <- scorepair$score

  scorepair$points <- factor(typenams[scorepair$score],levels = unique(typenams))
  #Levels: site taxon centroid


  ####

  scorepair$weight[is.na(scorepair$weight)] <- 0.01 *size.centroids * min(scorepair$weight,na.rm= TRUE) # the centroids
  scorepair$sizeweight  <- scorepair$weight

  #todo
  scorepair$sizeweight[scorepair$score%in% c("centroids","centroids_traits")] <-0.5 * stats::median(scorepair$sizeweight,na.rm= TRUE)
  scorepair$smoothweight  <-scorepair$weight
  scorepair$smoothweight[scorepair$score%in% c("centroids","centroids_traits")] <- NA# 0.0000001*min(scorepair$smoothweight,na.rm = TRUE)

  # missing for dcCA axis for centroids when smoothing
  scorepair$xforsmooth <- scorepair$dcCA1
  scorepair$xforsmooth[scorepair$score%in% c("centroids","centroids_traits")] <-   NA


  if(length(envfactor1)>1) {
    envlevels <- envlevels1 <- levels(envfactor1); #envfactor1 <- envfactor
    envlevels2 <- NULL
    if (sum(scorepair$score%in% c("centroids"))==0) envlevels0 <- NULL else envlevels0 <- envlevels
  } else {
    envlevels <-   NULL
    envlevels1 <- envlevels2 <-  "sites";
    if (sum(scorepair$score%in% c("centroids"))==0)
      envlevels0 <- NULL else envlevels0 <-scorepair$label[scorepair$score%in% c("centroids")]
    envfactor1 <- factor(rep("sites", nrow(object$data$dataEnv)))
  }


  if (length(traitfactor1)>1) {
    traitlevels <- traitlevels1 <- levels(traitfactor1); #traitfactor1 <- traitfactor;
    traitlevels2 <-  NULL
    if (sum(scorepair$score%in% c("centroids_traits"))==0) traitlevels0 <- NULL else traitlevels0 <- traitlevels
  }else {
    traitlevels <- NULL
    traitlevels1 <-   traitlevels2 <- "species";
    if (sum(scorepair$score%in% c("centroids_traits"))==0) traitlevels0 <- NULL else traitlevels0 <-scorepair$label[scorepair$score%in% c("centroids_traits")]
    traitfactor1 <- factor(rep("species", nrow(object$data$dataTraits)))
  }


  scorepair$groups  <- factor(c(envlevels1[envfactor1], envlevels0,traitlevels1[traitfactor1], traitlevels0),

                              levels = c(envlevels,traitlevels, envlevels2 ,traitlevels2))
  scorepair$names  <- rownames(scorepair)

  scorepair$centroidnames  <- scorepair$label
  scorepair$centroidnames[ ! scorepair$score%in%c("centroids","centroids_traits")] <- ""

  scorepair$centroidnames[scorepair$score%in%c("centroids_traits")] <-
        newname.list$centroidnames[[1]]
  scorepair$centroidnames[scorepair$score%in%c("centroids")] <-
    newname.list$centroidnames[[2]]

  names(scorepair)[1] <- paste("dcCA", axis ,sep = "")
  attr(scorepair, "condition") <- c(traitINcondition=traitINcondition,envINcondition=envINcondition)
  attr(scorepair, "levels") <- list(traitlevels=traitlevels,envlevels=envlevels)


res <- list(CWM_SNC = scorepair, trait_env_scores = trait_env_scores, newname.list = newname.list )
return(res)
}

setnames <- function(mod_scores, newnames = NULL){
  # mod_scores tidy = TRUE
  idTFt <- mod_scores$score =="intra_set_correlation_traits"
  idTFe <- mod_scores$score =="intra_set_correlation"
  oldnames <- list(traits = mod_scores$label[idTFt], env = mod_scores$label[idTFe])
  if (is.null(newnames)) {
    newnames <- oldnames
  } else {
    if (!length(newnames[[1]])==sum(idTFt)){
      print(paste("current trait names are:", paste0( mod_scores$label[idTFt], collapse = ",")))
      if (!length(newnames[[2]])==sum(idTFe)) {
      print(paste("current environmental names are:", paste0( mod_scores$label[idTFe], collapse = ",")))
      print(paste("newnames[[2]] should have length ", sum(idTFe), "; current length is", length(newnames[[2]])), sep ="")
      }
     # print(paste0("current trait names are:", mod_scores$label[idTFt], collapse = "," ))
     # stop(paste("newnames[[1]] should have length", length(idTFt)))
      stop(paste("newnames[[1]] should have length ", sum(idTFt), "current length is", length(newnames[[1]])),sep ="")

    }
    if (!length(newnames[[2]])==sum(idTFe)) {
      print(paste("current environmental names are:", paste0( mod_scores$label[idTFe], collapse = ",")))
      stop(paste("newnames[[2]] should have length ", sum(idTFe), "; current length is", length(newnames[[2]])), sep ="")
    }
  }
  idTFt <- mod_scores$score =="centroids_traits"
  idTFe <- mod_scores$score =="centroids"

  oldcentroidnames <- list(traits = mod_scores$label[idTFt], env = mod_scores$label[idTFe])
  # centroid names
  centroidnames <- list()
  ids <- which(oldnames[[1]] %in% oldcentroidnames[[1]])
  centroidnames[["traits"]] <- newnames[[1]][ids]
  ids <- which(oldnames[[2]] %in% oldcentroidnames[[2]])
  centroidnames[["env"]] <- newnames[[2]][ids]
  # regression names
  regnames <- list()
  nam_regr <- mod_scores$label[mod_scores$score == c("regression_traits")]
  # the variable name that is not in nam regr
  nam_not_regr <- oldnames[[1]][which(  !oldnames[[1]]%in% nam_regr)]# often the first of nam_centroids ...
  #name_not_among the regression coefs.
  idn <- oldnames[[1]][oldnames[[1]] %in% nam_not_regr]
  regnames[["traits"]] <- newnames[[1]][!oldnames[[1]]%in% idn]

  # reg env names  ---------------------------------------------------------------
  nam_regr <- mod_scores$label[mod_scores$score == c("regression")]
  # the variable name that is not in nam regr
  nam_not_regr <- oldnames[[2]][which(  !oldnames[[2]]%in% nam_regr)]# often the first of nam_centroids ...
  #name_not_among the regression coefs.
  idn <- oldnames[[2]][oldnames[[2]] %in% nam_not_regr]
  regnames[["env"]] <- newnames[[2]][!oldnames[[2]]%in% idn]
  return(newnames = list(newnames = newnames, weightnames = regnames, centroidnames = centroidnames))
}
