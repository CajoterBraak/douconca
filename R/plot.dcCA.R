#' Plot a single dc-CA axis with CWMs, SNCs, trait and environment scores.
#' 
#' @description
#' \code{plot.dcca} plots the CWMs and SNCs of a dc-CA axis against this axis,
#' with optional centroids and colors for groups of sites and/or species if 
#' available in the data.
#' 
#' @inheritParams getPlotdata
#' 
#' @param ... unused.
#' @param gradient_description character or 2-character vector for the trait
#' and environmental gradient, respectively specifying what to plot in the 
#' vertical line plots to describe the dc-CA axis (trait and environmental 
#' gradients). Default: \code{correlation} for intra-set correlations of both 
#' sets of variables with their dc-CA axis. Other values are: 
#' \code{c("regression", "weights", "tvalues", "inter_set_correlation")} for 
#' regression weights, t-values and (other) namely inter-set correlation, 
#' being the correlation of the SNCs and CWMs with the traits and environmental
#' variables, respectively.
#' @param facet logical. Default \code{TRUE} for CWMs and SNCs plots in 
#' separate panels. This parameter changes the position of the centroid 
#' names (from left to right for the environmental centroids). If 
#' \code{facet = FALSE} and \code{with_lines = TRUE}, the line fits ignore 
#' groups of species and of sites.
#' @param with_lines logical. Default \code{TRUE} for straight lines through
#' groups of points.
#' @param nspecies integer. Default \code{20} for including a vertical species 
#' plot with at most \code{nspecies} that have the highest contribution.
#' @param species_groups name of a variable in \code{dataTraits} of 
#' \code{\link{dc_CA}}. Default \code{NULL} for no grouping. NOT yet 
#' implemented.
#' @param verbose logical. Default \code{TRUE} for plotting the result.
#' @param widths relative widths of the CWM-SNC plot, the correlation/weight
#' plot and the species plot. (see \code{\link[gridExtra]{grid.arrange}}). 
#' Default \code{c(5, 1, 1)}.
#' 
#' @details
#' If you want to set new names, look at the names with all arguments default, 
#' i.e. \code{myplot <- plot(x)}, and then consult 
#' \code{myplot$nameList$newnames} for the order of the names of traits and
#' environmental variables. Note that covariates should not be in the list of
#' names. Contribution (in the definition of species selection in 
#' \code{nspecies}) is defined (as in CA) as the total species abundance in 
#' the (possibly, closed) data multiplied by the square of the score on 
#' the axis.
#'
#' If the \code{plot.dcca} returns the error \code{"Error in grid.Call"}, 
#' enlarge the plotting area or use \code{verbose = FALSE} and assign the 
#' result.
#'
#' @example demo/dune_plot_dcCA.R
#' @export
plot.dcca <- function(x, 
                      ...,
                      axis = 1,
                      gradient_description = "correlation",
                      envfactor = NULL, 
                      traitfactor = NULL, 
                      nspecies = 20, 
                      species_groups = NULL, 
                      widths = c(5, 1, 1),
                      newnames = NULL, 
                      facet = TRUE, 
                      remove_centroids = FALSE, 
                      with_lines = TRUE, 
                      verbose = TRUE) {
  if (!inherits(x, "dcca")) {
    stop("x should be of class dcca.\n")
  }
  stats_vals = c("regression", "weights", "correlations", "tvalues", 
                 "inter_set_correlation")
  if (length(gradient_description) == 1) {
    gradient_description <- match.arg(gradient_description[1],
                                      choices = stats_vals)
    gradient_description <- c(gradient_description, gradient_description)
  } else{
    gradient_description[1] <- match.arg(gradient_description[1],
                                         choices = stats_vals)
    gradient_description[2] <- match.arg(gradient_description[2],
                                         choices = stats_vals)
  }
  if (nspecies == 0) {
    widths <- c(widths[1], sum(widths[-1]))
  }
  pd <- getPlotdata(x, axis = axis, envfactor = envfactor, 
                    traitfactor = traitfactor, facet = facet, 
                    newnames = newnames, remove_centroids = remove_centroids)
  CWM_SNC <- plot_dcCA_CWM_SNC(x, axis = axis, envfactor = envfactor, 
                               traitfactor = traitfactor, facet = facet,
                               remove_centroids = remove_centroids, 
                               with_lines = with_lines, 
                               getPlotdata2plotdCCA = pd)
  trait_env_scores <- pd$trait_env_scores
  trait_env_scores$score <- factor(trait_env_scores$score)
  # trait score
  stats_scores <- list(traits = gradient_description[1], 
                       env = gradient_description[2])
  stats_scores <- lapply(stats_scores, function(x) {
    if (x == "regression") x <- "weights"
    return(x)
  })
  ylab_traits <- "composite trait"
  newnames <- "newnames"
  ncovariates <- 0
  if (stats_scores[[1]][1] == "weights") {
    idTF <- pd$trait_env_scores$score == "regression_traits"
    ncovariates <- sum(idTF) - length(pd$newNameList$weightnames$traits)
    if (ncovariates < 0) {
      print(pd$newNameList$weightnames$traits)
    }
    trait_title <- "Weight"
    newnames <- "weightnames"
  } else if (stats_scores[[1]][1] == "tvalues") {
    idTF <- pd$trait_env_scores$score == "t_values_traits"
    ncovariates <- sum(idTF) - length(pd$newNameList$weightnames$traits)
    if (ncovariates < 0) {
      print(pd$newNameList$weightnames$traits)
    }
    trait_title <- "t-value"
    newnames <- "weightnames"
  } else if (stats_scores[[1]][1] == "correlations") { # intra-set correlation
    idTF <- pd$trait_env_scores$score == "intra_set_correlation_traits"
    trait_title <- "correlation"
  } else { # inter-set correlation
    idTF <- pd$trait_env_scores$score == "correlation_traits"
    trait_title <- "corr. with SNC"
    ylab_traits <- "with SNC of dc-CA axis"
  }
  trait_scores <- pd$trait_env_scores[idTF, ]
  if (ncovariates > 0) {
    trait_scores <- trait_scores[-seq_len(ncovariates), , drop = FALSE]
  }
  trait_scores$label <- pd$newNameList[[newnames]]$traits
  ylab_env <-  "dc-CA axis"
  newnames <- "newnames"
  ncovariates <- 0
  if (stats_scores[[2]][1] == "weights") {
    idTF <- pd$trait_env_scores$score == "regression"
    ncovariates <- sum(idTF) - length(pd$newNameList$weightnames$env)
    if (ncovariates < 0) {
      print(pd$newNameList$weightnames$env)
    }
    env_title <- "weight"
    newnames <- "weightnames"
  } else if (stats_scores[[2]][1] == "tvalues") {
    idTF <- pd$trait_env_scores$score == "t_values"
    ncovariates <- sum(idTF) - length(pd$newNameList$weightnames$env)
    if (ncovariates < 0) {
      print(pd$newNameList$weightnames$traits)
    }
    env_title <- "t-value"
    newnames <- "weightnames"
    #ylab_env <-  "t-value in dc-CA axis"
  } else if (stats_scores[[2]][1] == "correlations") { # intra-set correlation
    idTF <- pd$trait_env_scores$score == "intra_set_correlation"
    env_title <- "correlation"
  } else { # inter-set correlation
    idTF <- pd$trait_env_scores$score == "correlation"
    env_title <- "corr. with CWM"
    ylab_env <- "with CWM of composite trait)"
  }
  env_scores <- pd$trait_env_scores[idTF, ]
  if (ncovariates > 0) {
    env_scores <- env_scores[-seq_len(ncovariates), , drop = FALSE]
  }
  env_scores$label <- pd$newNameList[[newnames]]$env
  namaxis <- names(env_scores)[1]
  if (!stats_scores[[1]][1] %in% c("weights", "tvalues", "correlations")) {
    trait_title <- "correlation"
  }
  if (gradient_description[1] == gradient_description[2]) {
    env_title <- ""
  } else {
    if (!stats_scores[[2]][1] %in% c("weights", "tvalues", "correlations")) {
      env_title <- "correlation"
    }
  }
  if (stats_scores[[1]][1] == "tvalues") {
    y_lab_interval <- 1 
  } else { 
    y_lab_interval <- 0.2
  }
  plot_traits <- plot_species_scores_bk(
    species_scores = trait_scores,
    ylab = ylab_traits,
    threshold = 0,
    y_lab_interval = y_lab_interval,
    speciesname = "label",
    scoresname = namaxis,
    selectname = "Fratio1",
    verbose = FALSE) + 
    ggplot2::ggtitle(trait_title)
  plot_env <- plot_species_scores_bk(
    species_scores = env_scores,
    ylab = ylab_env,
    threshold = 0,
    y_lab_interval = y_lab_interval,
    speciesname = "label",
    scoresname = namaxis,
    selectname = "Fratio1",
    verbose = FALSE) + 
    ggplot2::ggtitle(env_title)
  # species vertical plot
  plot_species <- fplot_species(pd, x, nspecies = nspecies, 
                                species_groups = species_groups)
  # plot arrange
  if (nspecies) {
    layout <- rbind(c(1, 2, 4), c(1, 3, 4))
    gg_object <- gridExtra::arrangeGrob(CWM_SNC, plot_traits, plot_env, 
                                        plot_species, layout_matrix = layout, 
                                        widths = widths)
  } else {
    layout<- rbind(c(1, 2), c(1, 3))
    gg_object <- gridExtra::arrangeGrob(CWM_SNC, plot_traits, plot_env, 
                                        layout_matrix = layout, 
                                        widths = widths)
  }
  # plot
  if (verbose) {
    tt <- try(suppressWarnings(gridExtra::grid.arrange(gg_object)))
    if (inherits(tt, "try-error")) {
      warning("Enlarge the plot area.\n")
    }
  }
  out <- list(plot = gg_object, nameList = pd$newNameList,
              separateplots = list(CWM_SNC = CWM_SNC, traits = plot_traits, 
                                   env = plot_env, species = plot_species))
  invisible(out)
}


#' @noRd
#' @keywords internal
fplot_species <- function(pd,
                          x, 
                          nspecies = 0, 
                          species_groups = NULL) {
  if (nspecies) {
    composite_trait <- pd$CWM_SNC[pd$CWM_SNC$score == "constraints_species", 1]
    contribution <- x$weights$columns * composite_trait ^ 2
    # just for later add a grouping
    if (!is.null(species_groups)) {
      # if one name in dataTraits take it.
      if (length(species_groups) == 1) {
        if (species_groups %in% names(x$data$dataTraits)) {
          species_groups <- x$data$dataTraits[[species_groups]] 
        } else {
          warning("species_groups not in names of dataTraits; ", 
                  "no grouping in the plot.\n")
          species_groups <- NULL
        }
      }
    }
    SNC_LC_mat <- cbind(composite_trait, contribution, species_groups)
    colnames(SNC_LC_mat) <- 
      c("composite_trait", "contribution", "species_group")[seq_len(ncol(SNC_LC_mat))]
    rownames(SNC_LC_mat) <- colnames(x$data$Y)
    sspecies <- sort(SNC_LC_mat[, "contribution"], decreasing = TRUE)
    threshold <- SNC_LC_mat[names(sspecies)[nspecies + 1], "contribution"]
    # ready for plotting
    plot_species <- plot_species_scores_bk(
      species_scores = SNC_LC_mat,
      ylab = "trait composite",
      threshold = threshold,
      y_lab_interval = 0.5,
      expand = 0.1,
      speciesname = NULL,
      scoresname = "composite_trait",
      selectname = "contribution",
      verbose = FALSE
    )
  } else {
    plot_species <- NULL
  }
  return(plot_species)
}
