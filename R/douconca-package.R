#' @title The package douconca performs double constrained correspondence 
#' analysis for trait-environment analysis in ecology
#'
#' @description
#' Double constrained correspondence analysis (dc-CA) for analyzing 
#' (multi-)trait (multi-)environment ecological data using library \code{vegan}
#' and native R code. It has a \code{formula} interface which allows to assess, 
#' for example, the importance of trait interactions in shaping ecological 
#' communities. The function \code{dc_CA} has an option to divide the abundance
#' data of a site by the site total, giving equal site weights. This division 
#' has the advantage that the multivariate analysis corresponds with an 
#' unweighted (multi-trait) community-level analysis, instead of being weighted.
#'
#' Throughout the two step algorithm of ter Braak et al. (2018) is used. This 
#' algorithm combines and extends community- (sample-) and species-level 
#' analyses, \emph{i.e.} (1) the usual community weighted means (CWM) 
#' regression analysis and (2) the species-level analysis of species-niche 
#' centroids (SNC) regression analysis. The SNC is the center of the realized 
#' niche of the species along an environmental variable or, in the case of 
#' dc-CA, an environmental gradient, \emph{i.e.} the dc-CA ordination axis.
#' Computationally, dc-CA can be carried out by a single singular value 
#' decomposition (ter Braak et al. 2018), but it is here computed in two steps.
#'
#' The first step uses canonical correspondence analysis 
#' (\code{\link[vegan]{cca}}) to regress the (transposed) abundance data on to 
#' the traits and the second step uses weighed redundancy analysis 
#' (\code{\link{wrda}} or, with equal site weights, \code{\link[vegan]{rda}})
#' to regress the CWMs of the orthonormalized traits, 
#' obtained from the first step, on to the environmental predictors. 
#' The second step is thus a community-level analysis.
#'
#' If \code{divideBySiteTotals = FALSE}, the second step uses 
#' \code{\link{wrda}} and performs a weighted redundancy analysis of the CWMs 
#' on to the environmental variables.
#'
#' Division of the abundance data by the site totals has the advantage that 
#' the resulting analysis (without dimension reduction, \emph{i.e.} retaining 
#' all dc-CA axes) corresponds with a series of unweighted community-level 
#' analyses, instead of the analyses being weighted.
#'
#' Warning: The \code{dcCA} package was built from \code{vegan} version 2.6-4 
#' and uses some of the internal structure of the \code{vegan} 
#' \code{\link[vegan]{cca.object}} in the not-exported functions 
#' \code{f_inertia} and \code{get_QR} in the source code file 
#' \code{functions_using_cca_object_internals.r}.
#'
#' The main user-functions are \code{\link{dc_CA}}, \code{\link{plot.dcca}}, 
#' \code{\link{scores.dcca}}, \code{\link{print.dcca}} and 
#' \code{\link{anova.dcca}}.
#'
#' @references
#' ter Braak, CJF, Šmilauer P, and Dray S. 2018. Algorithms and biplots for
#' double constrained correspondence analysis.
#' Environmental and Ecological Statistics, 25(2), 171-197.
#' \doi{10.1007/s10651-017-0395-x} 
#'
#' Oksanen, J., et al. (2022)
#' vegan: Community Ecology Package. R package version 2.6-4.
#' \url{https://CRAN.R-project.org/package=vegan}.
#'
#' @seealso \code{\link[vegan]{cca}} and \code{\link[vegan]{rda}}
#' 
#' @aliases douconca-package
#' @name douconca-package
#' @keywords internal
#' 
#' @importFrom stats as.formula contrasts delete.response median model.frame model.matrix reformulate terms
#' @importFrom rlang .data
"_PACKAGE"
NULL
