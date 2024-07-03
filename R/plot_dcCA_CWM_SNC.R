#' Plot the CWMs and SNCs of a single dc-CA axis.
#' 
#' @description
#' \code{plot_dcCA_CWM_SNC} plots the CWMs and SNCs of a dc-CA axis against 
#' this axis, with optional centroids and colors for groups of sites and/or 
#' species if available in the data.
#' 
#' @inheritParams getPlotdata
#' 
#' @param facet logical. Default \code{TRUE} for CWMs and SNCs plots in 
#' separate panels. This parameter changes the position of the centroid names
#' (from left to right for the environmental centroids). If \code{facet = TRUE}
#' and \code{with_lines = TRUE}, the line fits ignore groups of species and 
#' of sites.
#' @param with_lines logical. Default \code{TRUE} for straight lines through 
#' groups of points. \code{traitfactor = NA} and \code{envfactor = NA}. 
#' Centroids are not displayed in this case.
#' @param getPlotdata2plotdCCA the results of an \code{\link{getPlotdata}}. 
#' Default \code{NULL}.
#' 
#' @details
#' The argument \code{getPlotdata2plotdCCA} is to allow some modifications of
#' the data frame resulting from \code{\link{getPlotdata}}. The variable names
#' and score levels should remain untouched. \code{plot_dcCA_CWM_SNC} uses the
#' variables: \code{"dcCA\emph{k}"} with axis number \emph{k} and
#' \code{"CWM-SNC", "groups", "points", "sizeweight"} for the y-axis, coloring, 
#' shape and size of items, respectively.
#' 
#' @example demo/dune_plot_dcCA.R
#' 
#' @export
plot_dcCA_CWM_SNC <- function(object, 
                              axis = 1,
                              envfactor = NULL,
                              traitfactor = NULL,
                              facet = TRUE,
                              newnames = NULL,
                              remove.centroids = FALSE,
                              with_lines = TRUE,
                              getPlotdata2plotdCCA = NULL) {
  if (is.null(getPlotdata2plotdCCA)) {
    scorepair<- getPlotdata(object, axis = axis, envfactor = envfactor, 
                            traitfactor = traitfactor, newnames = newnames,
                            facet = facet, 
                            remove.centroids = remove.centroids)$CWM_SNC
  } else {
    scorepair = getPlotdata2plotdCCA$CWM_SNC
  }
  # plot
  namaxis <- paste0("dcCA", axis)
  p <- ggplot2::ggplot(data = scorepair, 
                       ggplot2::aes(x = .data[[namaxis]], 
                                    y = .data[["CWM-SNC"]],
                                    group = .data[["groups"]], 
                                    color = .data[["groups"]],
                                    shape = .data[["points"]], 
                                    size = .data[["sizeweight"]])) +
    ggrepel::geom_label_repel(ggplot2::aes(label = .data[["centroidnames"]]),
                              max.overlaps = 100) +
    ggplot2::geom_point() +
    ggplot2::ylab("") +
    ggplot2::xlab(paste("dc-CA axis", axis)) +
    ggplot2::scale_size(guide = "none")
  if (facet) {
    p <- p + ggplot2::facet_grid(type~., switch = "y" ) +
      ggplot2::scale_y_continuous(position = "right") 
  } else {
    p <- p + 
      ggplot2::ylab("CWM and SNC")
  }
  if (with_lines) {
    if (facet) {
      p <- p +
        ggplot2::geom_smooth(ggplot2::aes(x = .data[["xforsmooth"]],
                                          y = .data[["CWM-SNC"]],
                                          group = .data[["groups"]],
                                          weight = .data[["weight"]]),
                             linewidth = 1, method = "lm", formula = y ~ x,
                             na.rm = TRUE, inherit.aes = FALSE)
    } else {
      p <- p +
        ggplot2::geom_smooth(ggplot2::aes(x = .data[["xforsmooth"]], 
                                          y = .data[["CWM-SNC"]],
                                          group = .data[["type"]],
                                          weight = .data[["weight"]]),
                             linewidth = 1, method = "lm", formula = y ~ x,
                             na.rm = TRUE, inherit.aes = FALSE)
    }
  }
  TraitEnvINcondition <- attr(scorepair, "condition")
  TraitEnvLevels <- attr(scorepair, "levels")
  ltraits <- length(TraitEnvLevels[[1]])
  lenv <- length(TraitEnvLevels[[2]])
  if (ltraits > 0 && lenv > 0) {
    p <- p +
      ggplot2::guides(col = ggplot2::guide_legend(title = NULL, 
                                                  position = "bottom",
                                                  ncol = lenv, 
                                                  theme = ggplot2::theme(legend.byrow = TRUE)),
                      shape = ggplot2::guide_legend(title = NULL, 
                                                    position = "bottom", 
                                                    nrow = 2))
  } else if (ltraits == 0 && lenv == 0 && 
             sum(TraitEnvINcondition) == 0 && sum(remove.centroids) == 0) {
    p <- p + 
      ggplot2::scale_color_discrete(guide = "none") +
      ggplot2::guides(shape = ggplot2::guide_legend(title = NULL, 
                                                    position = "bottom", 
                                                    nrow = 1))
  } else {
    p <- p + 
      ggplot2::guides(col = ggplot2::guide_legend(title = NULL, 
                                                  position = "bottom", 
                                                  nrow = 1))
  }
  return(p)
}
