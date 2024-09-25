#' @title Vertical ggplot2 line plot of ordination scores
#'
#' @description
#' \code{plot_species_scores_bk} creates a vertical line plot of ordination 
#' scores with selection criterion for which scores to plot with names.
#' 
#' @param species_scores  a species-by-scores matrix, a data frame with 
#' row names (species names) or a tibble with variable with name
#' \code{speciesname} containing species names and a column or variable with 
#' name \code{scoresname} containing the scores (default: \code{"RDA1"}), e.g.
#' species scores from library \code{vegan}.
#' @param ylab y-axis label. Default: $b_k$.
#' @param y_lab_interval interval of the y-axis ticks. A tick at no effect (0) 
#' is always included; default: 0.5.
#' @param threshold species with criterion (specified by \code{selectname}) 
#' higher than the \code{threshold} are displayed. Default: 7 (which is the 
#' threshold F-ratio for testing a single regression coefficient at
#' \code{p = 0.01} with \code{60} df for the error in a multiple regression
#' of each single species onto the condition and the ordination axis). If 
#' \code{selectname} is not in \code{species_scores}, the \code{threshold} is 
#' divided by \code{14}, so that the default is 0.5.
#' @param speciesname name of the variable containing the species names 
#' (default \code{NULL} uses row names).
#' @param scoresname name of the column or variable containing the species 
#' scores to be plotted (default \code{"RDA1"}).
#' @param selectname name of the column or variable containing the criterion 
#' for the selection of species to be displayed. Default: \code{"Fratio1"}; 
#' if \code{selectname} is not found in \code{species_scores}, set 
#' to \code{scoresname}.
#' @param speciesgroup name of the factor, the levels of which
#' receive different colors in the vertical plot.
#' @param expand amount of extension of the line plot (default 0.2).
#' @param verbose logical for printing the number of species with names out of
#' the total number (default: \code{TRUE}).
#' 
#' @details
#' The absolute value of the criterion values is taken before selection, so 
#' that also the species scores of the ordination can serve as a criterion
#' (e.g. \code{selectname = "RDA1"}). The function has been copied from 
#' the \code{PRC} package at \url{https://github.com/CajoterBraak/PRC}.
#' 
#' The function is used in \code{\link{plot.dcca}}.
#'
#' @return a ggplot object
#' 
#' @example demo/dune_plot_gd.r
#' 
#' @export
plot_species_scores_bk <- function(species_scores, 
                                   ylab = "scores", 
                                   threshold = 7, 
                                   y_lab_interval = 0.5,
                                   speciesname = NULL, 
                                   scoresname = "RDA1",
                                   selectname = "Fratio1",
                                   speciesgroup = NULL,
                                   expand = 0.2, 
                                   verbose = TRUE) {
  # species_scores is a matrix or dataframe  with rownames and a column with 
  #  name scoresname (default: "RDA1") (species scores from vegan, for example)
  if (!is.null(speciesname) && speciesname %in% names(species_scores)) {
    sppnames <- species_scores[, speciesname]
  } else if (inherits(species_scores, c("matrix", "data.frame")) && 
             !is.null(rownames(species_scores))) {
    sppnames <- rownames(species_scores)
  } else {
    sppnames <- 1:nrow(species_scores)
  }
  if (is.matrix(species_scores)) {
    namcols <- colnames(species_scores) 
  } else { 
    namcols <- names(species_scores)
  }
  if (scoresname %in% namcols) {
    scores <- species_scores[, scoresname] 
  } else {
    stop("no score column found in plot_species_scores.\n")
  }
  if (selectname %in% namcols) {
    selectcrit <- species_scores[, selectname] 
  } else {
    selectcrit <- abs(scores)
    threshold <- threshold / 14
  }
  if (!is.null(speciesgroup) && speciesgroup %in% namcols) {
    species_groups <- factor(species_scores[, speciesgroup])
  } else {
    species_groups <- NA
  }
  species <- data.frame(species = sppnames, scores = scores,
                        selectcrit = selectcrit, species_groups = species_groups)
  speciesSub <- subset(species, abs(selectcrit) > threshold)
  speciesSub <- speciesSub[order(speciesSub$scores), ]
  ymin <- min(species$scores) * 1.1
  ymax <- max(species$scores) * 1.1
  fbreaks <- function(ymax, yInterval) {
    yInterval <- abs(yInterval)
    if (ymax > 0) {
      if (ymax > yInterval)  {
        br <- seq(from = yInterval, to = ymax, by = yInterval)
      } else {
        br <- yInterval
      }
    } else { # ymax <= 0
      if (ymax < -yInterval) {
        br <- rev(seq(from = -yInterval, to = ymax, by = -yInterval))
      } else {
        br <- -yInterval
      }
    }
    return(br)
  }
  if (ymin * ymax < 0) { #ymin <0, ymax >0
    bkBreaks <- round(c(fbreaks(ymin, y_lab_interval), 0, 
                        fbreaks(ymax, y_lab_interval)), 1)
  } else {
    if (ymin > 0) {
      bkBreaks <- round(c(0, fbreaks(ymax, y_lab_interval)), 1) 
    } else {
      bkBreaks <- round(c(fbreaks(ymin, y_lab_interval), 0), 1)
    }
  }
  speciesSub$y.label.loc <- seq(from = ymin, to = ymax,
                                length.out = nrow(speciesSub))
  labellineFull <- data.frame(species = rep(rownames(species), each = 2),
                              x.coor = rep(c(0.00, 0.02), nrow(species)),
                              y.coor = rep(species$scores, each = 2))
  p <- ggplot2::ggplot(data = species, 
                       ggplot2::aes(y = scores,
                                    group = .data[["species"]],
                                    x = 0)) +
    ggplot2::coord_cartesian(xlim = c(0, 1), 
                             ylim = c(min(ymin, min(bkBreaks)) - expand,
                                      max(ymax, bkBreaks) + expand)) +
    ggplot2::geom_line(data = labellineFull,
                       ggplot2::aes(x = .data[["x.coor"]],
                                    y = .data[["y.coor"]],
                                    group = .data[["species"]]), 
                       linewidth = 1.5 * 5 / 14,
                       colour = "#F51F63") +      
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_line(colour = "grey50"),
      axis.ticks.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_text(),
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank()
    ) +
    ggplot2::scale_y_continuous(expand = c(0, 0), breaks = bkBreaks) +
    ggplot2::ylab(ylab)
  if (nrow(speciesSub) == 0) {
    warning("After thresholding, no named species left to display in ", 
            "plot_species_scores_bk (perhaps from plotPRC).\n")
  } else {
    if (verbose) {
      cat( nrow(speciesSub), " species with names out of", 
           nrow(species), "species\n")
    }
    labelline <- data.frame(species = rep(speciesSub$species, each = 4),
                            x.coor = rep(c(0, 0.04, 0.14, 0.16),
                                         nrow(speciesSub)),
                            y.coor = rep(speciesSub$scores, each = 4))
    labelline$y.coor[seq(from = 3, to = nrow(labelline), by = 4)] <-
      speciesSub$y.label.loc
    labelline$y.coor[seq(from = 4, to = nrow(labelline), by = 4)] <-
      speciesSub$y.label.loc
    p <- p + ggplot2::geom_text(data = speciesSub,
                                ggplot2::aes(label = .data[["species"]],
                                             x = 0.18, 
                                             y = .data[["y.label.loc"]],
                                             color = .data[["species_groups"]]),
                                hjust = 0,
                                size = 7 * (5 / 14), fontface = "italic") +
      ggplot2::guides(color = "none") +
      ggplot2::geom_line(data = labelline, 
                         ggplot2::aes(x = .data[["x.coor"]],
                                      y = .data[["y.coor"]],
                                      group = .data[["species"]]),
                         linewidth = 1 * 5 / 14, colour = "grey50")
  }
  return(p)
}
