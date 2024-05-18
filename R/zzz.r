#' @title The package douconca performs double constrained correspondence analysis for trait-environment analysis in ecology
#'
#' @description
#' The \code{douconca} package analyzes multi-trait multi-environment ecological data by
#' double constrained correspondence analysis (ter Braak et al. 2018) using \code{vegan} and native R code.
#' It has a \code{formula} interface which allows to assess, for example, the importance
#' of trait interactions in shaping ecological communities.
#' Throughout the two step algorithm of ter Braak et al. (2018) is used. This algorithm
#' combines and extends community- (sample-) and species-level analyses, \emph{i.e.}
#' the usual community weighted means (CWM)-based regression analysis and the
#' species-level analysis of species-niche centroids (SNC)-based regression analysis.
#' The SNC is the center of the realized niche of the species along an environmental variable or,
#' in the case of dc-CA, an environmental gradient, \emph{i.e.} the dc-CA ordination axis.
#' Computationally, dc-CA can be carried out by a single singular value decomposition (ter Braak et al. 2018),
#' but it is here computed in two steps.
#' The default \code{\link{dc_CA}} (\code{divide.by.site.totals==TRUE}) is \code{vegan}-based,
#' The first step uses \code{\link[vegan]{cca}} to regress the (transposed)
#' abundance data on to the traits
#' and the second step uses \code{\link[vegan]{rda}} to regress the
#' CWMs of the orthonormalized traits, obtained from the first step, on to the environmental predictors.
#' The second step is thus a community-level analysis.
#' The abundance data are divided by the sample total
#' (i.e. 'closed') in the vegan-based version. This
#' has the advantage that this multivariate analysis corresponds with an unweighted (multi-trait)
#' community-level analysis, instead of being weighted.
#' If \code{divide.by.site.totals==FALSE}, the second step uses \code{\link{wrda}} and performs a
#' weighted redundancy analyis of the CWMs on to the environmental variables.
#'
#' In \code{\link{anova_species}}, the first step extracts the SNCs with respect to the ordination axes
#' from an existing dc-CA analysis (the SNCs, more precisely, a rotated version thereof, could have been obtained from a \code{\link[vegan]{cca}}
#' to regress the abundance data on to the environmental variables), and then
#' applies a weighted redundancy analysis of these SNCs with the traits as predictors. The second step is thus
#' a species-level analysis, but the final result (eigenvalues/ordination axes) is identical with that of the
#' analyses steps in \code{\link{dc_CA}} in which the second step is a community-level analysis.
#' The vegan-based analysis is computationally efficient for community-level (site-based) permutation tests.
#' but the a factor 20 or so slower for species-level permutation tests.
#' The technical reason for the closure, is that \code{vegan} \code{\link[vegan]{rda}} cannot do a weighted analysis.
#'
#' Warning: The \code{dcCA} package was built from \code{vegan} version 2.6-4 and uses some of the
#' internal structure of the \code{vegan} \code{\link[vegan]{cca.object}}
#' in the not-exported functions \code{f_inertia} and \code{get_QR}
#' in the source code file \code{functions_using_cca_object_internals.r}.
#'
#' The main user-functions are \code{\link{dc_CA}}, \code{\link{plot_dcCA}}, \code{\link{scores.dcca}},
#'  \code{\link{print.dcca}} and \code{\link{anova.dcca}}.
#'
#'
#'
#' @references
#'
#' ter Braak, CJF, Šmilauer P, and Dray S. 2018. Algorithms and biplots for
#' double constrained correspondence analysis.
#' Environmental and Ecological Statistics, 25(2), 171-197.
#' https://doi.org/10.1007/s10651-017-0395-x or
#' http://rdcu.be/ETPh
#'
#' Oksanen, J., et al. (2022)
#' vegan: Community Ecology Package. R package version 2.6-4.
#' http://CRAN.R-project.org/package=vegan.
#'
#' @seealso \code{\link[vegan]{cca}} and \code{\link[vegan]{rda}}
#' @name douconca
NULL
.onLoad <- function(libname = find.package("douconca"), pkgname = "douconca"){
  # CRAN Note avoidance
  if (getRversion() >= "2.15.1")
    utils::globalVariables(  c(".data"))
  invisible()
}