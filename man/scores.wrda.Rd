% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scores.wrda.R
\name{scores.wrda}
\alias{scores.wrda}
\title{Extract results of a weighted reduncancy analysis (wrda)}
\usage{
\method{scores}{wrda}(
  x,
  choices = c(1, 2),
  display = c("all"),
  scaling = "sym",
  which_cor = "in model",
  tidy = FALSE,
  ...
)
}
\arguments{
\item{x}{object of class \code{"dcca"}, \emph{i.e.} result of \code{\link{dc_CA}}.}

\item{choices}{integer vector of which axes to obtain. Default: all dc-CA axes.}

\item{display}{a character vector, one or more of
\code{c("all","species","sites","sp", "wa", "lc","bp", "cor", "ic", "reg","tval", "cn")}.
The most items are as in \code{\link[vegan]{scores.cca}}, except \code{"cor"} and \code{"ic"}, for
inter-set and intra-set correlations, respectively, and \code{"tval"} for
the (over-optimistic) t-values of the regression coefficients.}

\item{scaling}{numeric (1,2 or 3) or character \code{"sites", "species" or "symmetric"}. Default: "sym".
Either site- (1) or species- (2) related scores are scaled by eigenvalues,
and the other set of scores is left unscaled,
or with 3 both are scaled symmetrically by square root of eigenvalues. Negative values are treated as the
corresponding positive ones by \code{abs(scaling)}. See also}

\item{which_cor}{character vector environmental variables names
in the data frames for which inter-set correlations must calculated.
Default: a character ("in_model") for all traits and variables in the model,
including collinear variables and levels.}

\item{tidy}{Return scores that are compatible with \code{ggplot2}:
all scores are in a single data.frame, score type is identified by factor variable \code{score},
the names by variable \code{label}, and species weights (in \code{\link{dc_CA}} are in variable \code{weight}.
See \code{\link[vegan]{scores.cca}}.}

\item{...}{Other arguments passed to the function (currently ignored).}
}
\value{
A data frame if \code{tidy = TRUE}, a matrix if a single item is asked for and a named list of matrices if more than one item
is asked for. The following names can be included: \code{c("sites",
"constraints_sites", "centroids", "regression", "t_values" ,"correlation","intra_set_correlation", "biplot",
"species")}. Each matrix has an attribute \code{"meaning"} explaining its meaning.
With \code{tidy=TRUE}, the resulting data frame has attributes \code{"scaling"} and \code{"meaning"}; the latter has
two columns: (1) name of score type and (2) its meaning, usage and interpretation.

An example of the meaning of scores in scaling \code{"symmetric"} with \code{display ="all"}:
\describe{
\item{ sites }{ CMWs of the trait axes (constraints species) in scaling 'symmetric' optimal for biplots and,
 almost so, for inter-site distances. }
 \item{ constraints_sites }{ linear combination of the environmental predictors and the covariates
  (making the ordination axes orthogonal to the covariates)'
   in scaling 'symmetric' optimal for biplots and, almost so, for inter-site distances. }
 \item{ regression }{ mean, sd, VIF, standardized regression coefficients and their
 optimistic t-ratio in scaling 'symmetric'. }
 \item{ t_values }{ t-values of the coefficients of the regression of the CWMs of
 the trait composite on to the environmental variables }
 \item{ correlation }{ inter set correlation,
  correlation between environmental variables and the sites scores (CWMs) }
 \item{ intra_set_correlation }{ intra set correlation, correlation between
  environmental variables and the dc-ca axis (constrained sites scores) }
 \item{ biplot }{ biplot scores of environmental variables for display with biplot-traits
 for fourth-corner correlations in scaling 'symmetric'. }
  \item{ centroids }{ environmental category means of the site scores in scaling 'symmetric'
 optimal for biplots and, almost so, for inter-environmental category distances.}
 \item{ species }{ SNC on the environmental axes (constraints sites) in scaling 'symmetric'
  optimal for biplots and, almost so, for inter-species distances. }
 \item{ constraints_species }{ linear combination of the traits and the trait covariates
  (making the ordination axes orthogonal to the covariates)' in scaling 'symmetric' optimal
   for biplots and, almost so, for inter-species distances. }
 \item{ regression_traits }{ mean, sd, VIF, standardized regression coefficients and
  their optimistic t-ratio in scaling 'symmetric'. }
 \item{ t_values_traits }{ t-values of the coefficients of the regression of the
  SNCs along a dc-CA axis on to the traits }
 \item{ correlation_traits }{ inter set correlation,
   correlation between traits and the species scores (SNCs) }
 \item{ intra_set_correlation_traits }{ intra set correlation,
  correlation between traits and the dc-ca axis (constrained species scores) }
 \item{ biplot_traits }{ biplot scores of traits for display with biplot scores
  for fourth-corner correlation in scaling 'symmetric'. }
  \item{ centroids_traits }{ trait category means of the species scores in scaling
   'symmetric' optimal for biplots and, almost so, for inter-trait category distances. }
}

The statements on optimality for distance interpretations are based on the \code{scaling}
and the relative magnitude of the dc-CA eigenvalues of the chosen axes.
}
\description{
This function works very much like the \code{vegan} \code{\link[vegan]{scores}} function,
in particular \code{\link[vegan]{scores.cca}}, with the additional results such
as regression coefficients and linear combinations of traits \code{('regr_traits','lc_traits')}
All scores from CA obey the so called transition formulas and so do the scores of CCA and dc-CA.
The differences are, for CCA, that the linear combinations of environmental variables
(the \emph{constrained} site scores)
replace the usual (\emph{unconstrained}) site scores, and for dc-CA,
that the linear combinations of traits (the \emph{constrained} species scores)
also replace the usual (\emph{unconstrained}) species scores.
}
\details{
The function is modeled after \code{\link[vegan]{scores.cca}}.

 If you get the error message:  'arg' should be one of "sites", "species", "both",
 then the vegan scores function has been called, instead of the one of douconca.
 The work-around is to use douconca::scores() instead of scores() only.

 An example of which_cor is: \code{which_cor = c("acidity","humidity")}
}
\examples{
data("dune_trait_env")
# rownames are carried forward in results
rownames(dune_trait_env$comm) <- dune_trait_env$comm$Sites
response = dune_trait_env$comm[, -1]  # must delete "Sites"
w <- rep(1,20); w[1:10]<- 8; w[17:20] <- 0.5

object <- wrda(formula= ~A1+Moist+Mag+Use+Condition(Manure),
          response = response, data= dune_trait_env$envir, weights = w)
print(object) # Proportions equal to those Canoco 5.15
mod_scores <- douconca::scores(object, display = "all")
douconca::scores(object, which_cor = c("A1", "Manure"), display = c("cor"))
anova(object)

# The default is equal weights, which allows checking against vegan
#object <- wrda(formula= ~A1+Moist+Mag+Use+Manure,
object1 <- wrda(formula= ~A1+Moist+Mag+Use+Condition(Manure),
               response = response, data= dune_trait_env$envir)
print(object1)

# compare with vegan::rda
#object2 <- vegan::rda(formula= response~A1+Moist+Mag+Use+Manure, data= dune_trait_env$envir)
object2 <- vegan::rda(formula= response~A1+Moist+Mag+Use+Condition(Manure), data= dune_trait_env$envir)
print(object2)
object1$CCA$eig/ object2$CCA$eig #
object1$tot.chi/object2$tot.chi
object1$CCA$tot.chi/object2$CCA$tot.chi
range(object1$CCA$u/ object2$CCA$u)

myconst <- sqrt(object1$Nobs*object2$tot.chi)
range(object1$site_axes$site_scores$lc_env_scores/ vegan::scores(object2, choices = seq_len(ncol(object1$site_axes$site_scores$lc_env_scores)),
                                                                display= "lc", scaling= "sites",const = myconst)) # 1 1

}
