#' Plot the CWMs and SNCs of a single dc-CA axis.
#' @description
#' \code{plot_dcCA_CWM_SNC} plots the CWMs and SNCs of a dc-CA axis against this axis,
#'  with optional centroids and colors for groups of sites and/or species if available in the data.
#' @inheritParams getPlotdata
#' @param facet logical. Default \code{TRUE} for CWMs and SNCs plots in separate panels.
#' This parameter changes the position of the centroid names (from left to right for the environmental
#' centroids).
#' If \code{facet = TRUE} and \code{with_lines = TRUE}, the line fits ignore groups of species and of sites.
#' @param with_lines logical. Default \code{TRUE} for straight lines though groups of points.
#' \code{traitfactor = NA} and \code{envfactor = NA}. Centroids are not displayed in this case.
#' @param getPlotdata2plotdCCA the results of an \code{\link{getPlotdata}}. Default \code{NULL}.
#' @details
#' The lines with \code{with_lines=TRUE} do no use the weights in this version and
#' give an extra band for an not-existing line (for missing centroids).
#' The argument \code{getPlotdata2plotdCCA} is to allow some modifications of the data frame resulting
#' from \code{\link{getPlotdata}}. The variable names and score levels should remain untouched.
#' \code{plot_dcCA_CWM_SNC} uses the variables: \code{"dcCA\emph{k}"} with axis number \emph{k} and
#' \code{"CWM-SNC", "groups", "points", "sizeweight"} for the y-axis, coloring, shape and size of
#' items, respectively.
#' @example demo/dune_plot_dcCA.R
#' @export
#'
plot_dcCA_CWM_SNC <- function(object, axis=1, envfactor=NULL, traitfactor=NULL, #size.centroids = 1,
         facet = TRUE, newnames = NULL, remove.centroids=FALSE, with_lines = TRUE, getPlotdata2plotdCCA=NULL){

# todo: outcomment if weight in smooth and extra band problems have been solved
if (facet) with_lines <- FALSE
if (is.null(getPlotdata2plotdCCA)){
scorepair<- getPlotdata(object, axis=axis, envfactor=envfactor, traitfactor=traitfactor,
           newnames = newnames,              # size.centroids=size.centroids,
           facet = facet, remove.centroids=remove.centroids)$CWM_SNC
} else scorepair = getPlotdata2plotdCCA$CWM_SNC
# plot ------------------------------------------------------

#cbbpalette <- c( "#D55E00","#303030", "#009E73", "#56B4E9",  "#F0E442", "#0072B2", "#CC79A7","black") # OK rev(blue, green, black)
namaxis <- paste("dcCA",axis, sep = "")
pp2<- ggplot2::ggplot(data= scorepair, ggplot2::aes(x = .data[[namaxis]], y = .data[["CWM-SNC"]],
                        group=.data[["groups"]], color = .data[["groups"]],
                        shape= .data[["points"]], size = .data[["sizeweight"]]))+
  ggrepel::geom_label_repel(ggplot2::aes(label = .data[["centroidnames"]]),max.overlaps = 100)+ggplot2::geom_point()+
  ggplot2::ylab("")+
  ggplot2::scale_size(guide="none")

if (facet) pp2 <- pp2+ ggplot2::facet_grid(type~., switch = "y" )+ggplot2::scale_y_continuous(position = "right")

#mylm <- function(fitweights, ...){stats::lm(formula= y~x, weights = fitweights,...)}

if (with_lines){
  #if (sum(remove.centroids)==0 && !facet)
  if (facet)  pp2 <-  pp2 +  ggplot2::geom_smooth(ggplot2::aes(x= .data[["xforsmooth"]]),
              linewidth = 1, method = stats::lm, na.rm=TRUE)#+ scale_shape(guide="none")# ,
                                        # method.args = list(weight = scorepair[["smoothweight"]]))
    else
    pp2 <-  pp2 +  ggplot2::geom_smooth(ggplot2::aes(x= .data[["xforsmooth"]],group=.data[["type"]]),linewidth = 1, method = stats::lm, na.rm=TRUE)
                           #  method.args = list(weight = scorepair[["smoothweight"]]))
}

TraitEnvINcondition <- attr(scorepair,"condition")
TraitEnvLevels <- attr(scorepair,"levels")
ltraits <-length(TraitEnvLevels[[1]])
lenv <-length(TraitEnvLevels[[2]])

if (ltraits >0 && lenv>0 ){
   pp3 <- pp2 +
     ggplot2::guides(col = ggplot2::guide_legend(title = NULL, position = "bottom",ncol=lenv, theme = ggplot2::theme(legend.byrow = TRUE) ),
            shape = ggplot2::guide_legend(title = NULL, position = "bottom", nrow = 2))#,
 } else if (ltraits==0 && lenv==0 && sum(TraitEnvINcondition)==0 && sum(remove.centroids)==0){
  pp3 <- pp2 + ggplot2::scale_color_discrete(guide="none")+
    ggplot2::guides(shape = ggplot2::guide_legend(title = NULL, position = "bottom", nrow = 1))
} else{
  pp3 <- pp2 + ggplot2::guides(col = ggplot2::guide_legend(title = NULL, position = "bottom", nrow = 1))
}
myxlab <- paste("dc-CA axis", axis)
return(suppressWarnings(pp3+ggplot2::xlab(myxlab)))
}
