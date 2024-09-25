#' Utility function: extracting data from a \code{\link{dc_CA}} object for 
#' plotting a single axis by your own code or \code{\link{plot.dcca}}.
#' 
#' @description
#' \code{getPlotdata} extracts data from a \code{\link{dc_CA}} object for 
#' plotting the CWMs and SNCs of a single axis.
#' 
#' @param x results from \code{\link{dc_CA}} of class \code{dcca}.
#' @param axis the axis number to get (default 1).
#' @param envfactor name of row factor to display as color and lines in the CWM
#' plot (default \code{NULL}). The default extracts the factor from the 
#' environmental model. If set to \code{NA}, no additional coloring and lines
#' are displayed in \code{\link{plot.dcca}}. The parameter sets the 
#' \code{groups} variable in the \code{CWM_SNC} data frame of the return 
#' value/in the plot.
#' @param traitfactor name of column factor to display as color and lines in
#' the SNC plot (default \code{NULL}). The default extracts the factor from
#' the trait model. If set to \code{NA}, no additional coloring and lines are
#' displayed in \code{\link{plot.dcca}}. The parameter sets the \code{groups} 
#' variable in the \code{CWM_SNC} data frame of the return value/in the plot.
#' @param newnames a list  with two elements: names for traits and for 
#' environmental variables, default \code{NULL} for names derived from the 
#' result of \code{\link{scores.dcca}} with \code{tidy = TRUE}.
#' @param remove_centroids logical to remove any centroids from the plot data 
#' (default \code{FALSE}). Can be a two-vector, \emph{e.g.} 
#' \code{c(TRUE, FALSE)} to remove only the environmental centroids.
#' @param facet logical. Default \code{TRUE} for CWMs and SNCs plots in 
#' separate panels. If \code{FALSE}, this parameter changes the position of 
#' the environmental centroid names (from left to right).
#' 
#' @return A list with three components
#' \describe{
#' \item{CWM_SNC}{a data.frame containing plot data}
#' \item{trait_env_scores}{a vector of scores per trait/environment}
#' \item{newNameList}{a vector of new names to be used in the plot}
#' }
#' 
#' @example demo/dune_plot_dcCA.R
#' 
#' @export
getPlotdata <- function(x, 
                        axis = 1,
                        envfactor = NULL,
                        traitfactor = NULL,
                        newnames = NULL,
                        facet = TRUE,
                        remove_centroids = FALSE) {
  size.centroids <- 1
  # getPlotdata function
  if (length(remove_centroids) == 1) {
    remove_centroids <- c(remove_centroids, remove_centroids)
  }
  traitINcondition <- envINcondition <- FALSE
  if (is.null(envfactor)) {
    ff <- get_Z_X_XZ_formula(x$formulaEnv, x$data$dataEnv)
    if (ff$formula_Z == ~1) {
      envfactor <- ff$focal_factor[1] 
      envINcondition <- FALSE 
    } else {
      envfactor <- ff$Condi_factor[1]
      envINcondition <- TRUE
    }
  }
  if (is.null(traitfactor)) {
    ff <- get_Z_X_XZ_formula(x$formulaTraits, x$data$dataTraits)
    if (ff$formula_Z == ~1) {
      traitfactor <- ff$focal_factor[1]
      traitINcondition <- FALSE
    } else {
      traitfactor <- ff$Condi_factor[1]
      traitINcondition <- TRUE
    }
  }
  if (is.null(traitfactor)) traitfactor <- NA
  if (is.null(envfactor)) envfactor <- NA
  # end of set env and traitfactor
  mod_scores <- scores(x, choices = axis, tidy = TRUE, 
                       scaling = "symmetric")
  if (!"species" %in% levels(mod_scores$score)) {
    # cannot do the species plot, can do CWM plot only
    stop("no unconstrained species scores, make a CWM plot instead.")
  }
  newNameList <- setnames(mod_scores, newnames = newnames)
  idTFc <- mod_scores$score %in% c("constraints_sites", "constraints_species", 
                                   "centroids", "centroids_traits")
  idTFu <- mod_scores$score %in% c("sites", "species", "centroids", 
                                   "centroids_traits")
  con_scores <- mod_scores[idTFc, ]
  uncon_scores <- mod_scores[idTFu, ]
  uncon_scores$score <- factor(uncon_scores$score)
  colnames(uncon_scores)[1] <- "CWM-SNC"
  trait_env_scores <- mod_scores[!(idTFc | idTFu), ]
  scorepair <- cbind(con_scores[, -4], uncon_scores[, c(1, 4), drop = FALSE])
  scorepair$score <- factor(scorepair$score)
  names(scorepair)[1] <- "dcCA1"
  # this removes centroids if traitfactor = NA
  if (remove_centroids[1]) {
    scorepair <- scorepair[scorepair$score != "centroids_traits", ]
  }
  if (remove_centroids[2]) {
    scorepair <-  scorepair[scorepair$score != "centroids", ]
  }
  if (!is.na(envfactor) || is.character(envfactor)) {
    if (is.na(envfactor)) { 
      envfactor1 <- envfactor 
    } else if(envfactor %in% names(x$data$dataEnv) ) {
      envfactor1 <- x$data$dataEnv[[envfactor]]
    } else {
      stop(envfactor, " must be in ", 
           paste0(names(x$data$dataEnv), collapse = ", "), ".\n")
    }
  } else {
    envfactor1 <- envfactor
  }
  if (!is.na(traitfactor) || is.character(traitfactor)) {
    if (is.na(traitfactor)) {
      traitfactor1 <- traitfactor
    } else if(traitfactor %in% names(x$data$dataTraits)) {
      traitfactor1 <- x$data$dataTraits[[traitfactor]]
    } else {
      stop(traitfactor, " must be in ", 
           paste0(names(x$data$dataTraits), collapse = ", "), ".\n")
    }
  } else {
    traitfactor1<- traitfactor
  }
  if (facet) {
    scorepair$dcCA1[scorepair$score%in% c("centroids", "centroids_traits")] <- 
      min(scorepair$dcCA1, na.rm = TRUE)
  } else {
    scorepair$dcCA1[scorepair$score == "centroids_traits"] <-
      min(scorepair$dcCA1, na.rm = TRUE)
    scorepair$dcCA1[scorepair$score == "centroids"] <-
      max(scorepair$dcCA1,na.rm = TRUE)
  }
  #### could be integrated in below....
  if (nlevels(scorepair$score) == 4) {
    typenams <- c("site", "taxon", "centroid", "centroid")
    scorepair$score <- factor(scorepair$score, 
                              levels = levels(scorepair$score)[c(3, 4, 1, 2)])
  } else if (nlevels(scorepair$score) > 2) {
    typenams <- c("site", "taxon", "centroid")
    scorepair$score <- factor(scorepair$score, 
                              levels = levels(scorepair$score)[c(2, 3, 1)])
  } else {
    typenams <- c("site", "taxon")
  }
  if (nlevels(scorepair$score) == 4) {
    typenams3 <- rep(c("CWM of trait composite", "SNC of dc-CA axis"), 2) 
  } else if (nlevels(scorepair$score) == 2) {
    typenams3 <- c("CWM of trait composite", "SNC of dc-CA axis") 
  } else if (nlevels(scorepair$score) == 3) {
    if ("centroids_traits" %in% levels(scorepair$score)) {
      typenams3 <- c("CWM of trait composite", rep("SNC of dc-CA axis", 2))
    } else {
      typenams3 <- c("CWM of trait composite", 
                     "SNC of dc-CA axis",
                     "CWM of trait composite")
    }
  }
  scorepair$type <-  typenams3[scorepair$score]
  scorepair$points <- factor(typenams[scorepair$score],
                             levels = unique(typenams))
  #Levels: site taxon centroid
  scorepair$weight[is.na(scorepair$weight)] <- 
    0.01 * size.centroids * min(scorepair$weight, na.rm = TRUE) # the centroids
  scorepair$sizeweight  <- scorepair$weight
  #todo
  scorepair$sizeweight[scorepair$score %in% c("centroids", "centroids_traits")] <-
    0.5 * median(scorepair$sizeweight,na.rm = TRUE)
  scorepair$smoothweight <- scorepair$weight
  scorepair$smoothweight[scorepair$score %in% c("centroids", "centroids_traits")] <- NA
  # missing for dcCA axis for centroids when smoothing
  scorepair$xforsmooth <- scorepair$dcCA1
  scorepair$xforsmooth[scorepair$score %in% c("centroids", "centroids_traits")] <- NA
  if (length(envfactor1) > 1) {
    envlevels <- envlevels1 <- levels(envfactor1)
    envlevels2 <- NULL
    if (sum(scorepair$score == "centroids") == 0) {
      envlevels0 <- NULL
    } else {
      envlevels0 <- envlevels
    }
  } else {
    envlevels <- NULL
    envlevels1 <- envlevels2 <- "sites"
    if (sum(scorepair$score == "centroids") == 0) {
      envlevels0 <- NULL 
    } else {
      envlevels0 <- scorepair$label[scorepair$score == "centroids"]
    }
    envfactor1 <- factor(rep("sites", nrow(x$data$dataEnv)))
  }
  if (length(traitfactor1) > 1) {
    traitlevels <- traitlevels1 <- levels(traitfactor1)
    traitlevels2 <- NULL
    if (sum(scorepair$score == "centroids_traits") == 0) {
      traitlevels0 <- NULL
    } else {
      traitlevels0 <- traitlevels
    }
  } else {
    traitlevels <- NULL
    traitlevels1 <- traitlevels2 <- "species"
    if (sum(scorepair$score == "centroids_traits") == 0) {
      traitlevels0 <- NULL 
    } else {
      traitlevels0 <- scorepair$label[scorepair$score == "centroids_traits"]
    }
    traitfactor1 <- factor(rep("species", nrow(x$data$dataTraits)))
  }
  scorepair$groups <- factor(c(envlevels1[envfactor1], 
                               envlevels0,traitlevels1[traitfactor1], 
                               traitlevels0),
                             levels = c(envlevels, traitlevels, 
                                        envlevels2, traitlevels2))
  scorepair$names  <- rownames(scorepair)
  scorepair$centroidnames <- scorepair$label
  scorepair$centroidnames[!scorepair$score %in% c("centroids", "centroids_traits")] <- ""
  scorepair$centroidnames[scorepair$score == "centroids_traits"] <-
    newNameList$centroidnames[[1]]
  scorepair$centroidnames[scorepair$score == "centroids"] <-
    newNameList$centroidnames[[2]]
  names(scorepair)[1] <- paste0("dcCA", axis)
  attr(scorepair, "condition") <- 
    c(traitINcondition = traitINcondition, envINcondition = envINcondition)
  attr(scorepair, "levels") <- 
    list(traitlevels = traitlevels, envlevels = envlevels)
  res <- list(CWM_SNC = scorepair, 
              trait_env_scores = trait_env_scores, 
              newNameList = newNameList)
  return(res)
}

#' @noRd
#' @keywords internal
setnames <- function(mod_scores, 
                     newnames = NULL){
  # mod_scores tidy = TRUE
  idTFt <- mod_scores$score == "intra_set_correlation_traits"
  idTFe <- mod_scores$score == "intra_set_correlation"
  oldnames <- list(traits = mod_scores$label[idTFt], 
                   env = mod_scores$label[idTFe])
  if (is.null(newnames)) {
    newnames <- oldnames
  } else {
    if (length(newnames[[1]]) != sum(idTFt)) {
      message("current trait names are:", 
              paste0(mod_scores$label[idTFt], collapse = ", "))
      if (length(newnames[[2]]) != sum(idTFe)) {
        message("current environmental names are:", 
                paste0( mod_scores$label[idTFe], collapse = ", "))
        warning("newnames[[2]] should have length ", sum(idTFe), 
                "; current length is", length(newnames[[2]]), ".\n")
      }
      stop("newnames[[1]] should have length ", sum(idTFt), 
           "current length is", length(newnames[[1]]), ".\n")
    }
    if (!length(newnames[[2]]) == sum(idTFe)) {
      message("current environmental names are:", 
              paste0( mod_scores$label[idTFe], collapse = ", "))
      stop("newnames[[2]] should have length ", sum(idTFe), 
           "; current length is", length(newnames[[2]]), ".\n")
    }
  }
  idTFt <- mod_scores$score == "centroids_traits"
  idTFe <- mod_scores$score == "centroids"
  oldcentroidnames <- list(traits = mod_scores$label[idTFt], 
                           env = mod_scores$label[idTFe])
  # centroid names
  centroidnames <- list()
  centroidnames[["traits"]] <- 
    newnames[[1]][oldnames[[1]] %in% oldcentroidnames[[1]]]
  centroidnames[["env"]] <- 
    newnames[[2]][oldnames[[2]] %in% oldcentroidnames[[2]]]
  # regression names
  regnames <- list()
  nam_regr <- mod_scores$label[mod_scores$score == "regression_traits"]
  # the variable name that is not in nam regr
  # often the first of nam_centroids ...
  nam_not_regr <- oldnames[[1]][!oldnames[[1]] %in% nam_regr]
  #name_not_among the regression coefs.
  idn <- oldnames[[1]][oldnames[[1]] %in% nam_not_regr]
  regnames[["traits"]] <- newnames[[1]][!oldnames[[1]] %in% idn]
  # reg env names  
  nam_regr <- mod_scores$label[mod_scores$score == "regression"]
  # the variable name that is not in nam regr
  # often the first of nam_centroids ...
  nam_not_regr <- oldnames[[2]][!oldnames[[2]] %in% nam_regr]
  #name_not_among the regression coefs.
  idn <- oldnames[[2]][oldnames[[2]] %in% nam_not_regr]
  regnames[["env"]] <- newnames[[2]][!oldnames[[2]] %in% idn]
  return(list(newnames = newnames, 
              weightnames = regnames,
              centroidnames = centroidnames))
}
