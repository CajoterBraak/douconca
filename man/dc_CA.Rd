% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dc_CA.r
\name{dc_CA}
\alias{dc_CA}
\title{Performs (weighted) double constrained correspondence analysis (dc-CA)}
\usage{
dc_CA(
  formulaEnv = NULL,
  formulaTraits = NULL,
  response = NULL,
  dataEnv = NULL,
  dataTraits = NULL,
  divide.by.site.totals = TRUE,
  dc_CA_object = NULL,
  verbose = TRUE
)
}
\arguments{
\item{formulaEnv}{formula or one-sided formula for the rows (samples) with row predictors in \code{dataEnv}.
When two-sided, the left hand side of the formula is not used.
Specify row covariates (if any ) by adding \code{+ Condition(covariate-formula)}
to \code{formulaEnv} as in \code{\link[vegan]{rda}}. The \code{covariate-formula} should not contain a \code{~} (tilde).
Default: \code{~.}, i.e. all variables in \code{dataEnv} are predictor variables.}

\item{formulaTraits}{formula or one-sided formula for the columns (species) with colum predictors in \code{dataTraits}.
When two-sided, the left hand side of the formula is not used. Specify column covariates (if any ) by  adding \code{+ Condition(covariate-formula)}
to \code{formulaTraits} as in \code{\link[vegan]{cca}}. The \code{covariate-formula} should not contain a \code{~} (tilde).
Default: \code{~.}, i.e. all variables in \code{dataTraits} are predictor traits.}

\item{response}{matrix, data frame of the abundance data (dimension \emph{n} x \emph{m}) or
object with community weighted means from \code{\link{fCWM_SNC}}. See Details for analyses starting from
community weighted means.
Rownames of \code{response}, if any, are carried through.}

\item{dataEnv}{matrix or data frame of the row predictors, with rows corresponding to those in \code{response}.
(dimension \emph{n} x \emph{p}).}

\item{dataTraits}{matrix or data frame of the column predictors,
with rows corresponding to the columns in \code{response}.(dimension \emph{m} x \emph{q}).}

\item{divide.by.site.totals}{logical;
default \code{TRUE} for closing the data by dividing the rows in the \code{response} by their total.}

\item{dc_CA_object}{optional object from an earlier run of this function. Useful if the
same formula for the columns (\code{formulaTraits}), \code{dataTraits} and \code{response} are used
with a new formula for the rows. If set, the data of the previous run is used and the result of its first step
is taken for the new analysis.
The \code{data$Y} element in the object can be set to \code{NULL}, which is useful,
In the case of non-public abundance data and a wish for reproducibility of the analysis
use the function \code{fCWM_SNC} to create an object for use as \code{response} argument
in a new call to \code{dc_CA}. From this object, the analysis can be reproduced from CWMs
(and, for full analysis, SNCs) instead of from abundance data.}

\item{verbose}{logical for printing a simple summary (default: TRUE)}
}
\value{
A list of \code{class} \code{dccav}; that is a list with elements
\describe{
\item{CCAonTraits}{a \code{\link[vegan]{cca.object}} from the \code{\link[vegan]{cca}} analysis
of the transpose of the closed \code{response} using formula \code{formulaTraits}. }
\item{formalaTraits}{the argument \code{formulaTraits}. If the
formula was \code{~.}, it was changed to explicit trait names.}
\item{data}{a list of \code{Y} (response data after removing empty rows and columns and after closure)
and \code{dataEnv} and \code{dataTraits}.}
\item{weights}{a list of unit-sum weights of row and columns.
The names of the list are \code{c("row","columns")}, in that order.}
\item{Nobs}{number of sites (rows).}
\item{CWMs_orthonormal_traits}{Community weighted means w.r.t. orthonormalized traits.}
\item{RDAonEnv}{a \code{\link[vegan]{cca.object}} from the \code{\link[vegan]{rda}} analysis
of the column scores of the \code{cca}, which are the CWMs of orthonormalized traits,
using formula \code{formulaEnv}. }
\item{formalaEnv}{the argument \code{formulaEnv}. If the
formula was \code{~.}, it was changed to explicit environmental variable names. }
\item{eigenvalues}{the dc-CA eigenvalues (same as those of the \code{\link[vegan]{rda}} analysis)}
\item{inertia}{a one-column matrix with four inertias (weighted variances):
\itemize{
\item total: the total inertia.
\item conditionT: the inertia explained by the condition in \code{formulaTraits}
if present (neglecting row constraints).
\item traits_explain: the inertia explained by the traits (neglecting the row predictors and any
condition in \code{formulaTraits}).
This is the maximum that the row predictors could explain in dc-CA
(the sum of the following two items is thus less than this value).
\item conditionE: the trait-constrained inertia explained by the condition in \code{formulaEnv}.
\item constraintsTE: the trait-constrained inertia explained by the predictors (without the row covariates).
}
 }
}
If \code{verbose} is \code{TRUE} (or after \code{out <-print(out)} is invoked )
there are three more items (in this version).
\itemize{
\item \code{c_traits_normed}: mean, sd, VIF and (regression) coefficients of
 the traits that define the dc-CA trait axes (composite traits), and their optimistic t-ratio.
\item \code{c_env_normed}:  mean, sd, VIF and (regression) coefficients of the environmental variables that define the dc-CA axes
 in terms of the environmental variables (composite gradients), and their optimistic t-ratio.
\item \code{species_axes}: a list with four items
 \itemize{
 \item \code{species_scores}: a list with names \code{c("species_scores_unconstrained", "lc_traits_scores")} with the
 matrix with species niche centroids along the dc-CA axes (composite gradients) and
 the matrix with linear combinations of traits.
 \item \code{correlation}: a matrix with inter-set correlations of the traits with their SNCs.
 \item \code{b_se}: a matrix with (unstandardized) regression coefficients for traits and their optimistic standard errors.
 \item \code{R2_traits}: a vector with coefficient of determination (R2) of the SNCs on to the traits.
 The square-root thereof could be called the species-trait correlation in analogy with
 the species-environment correlation in CCA.
 }
 \item \code{sites_axes}: a list with four items
 \itemize{
 \item \code{site_scores}: a list with names \code{c("site_scores_unconstrained", "lc_env_scores")} with the
 matrix with community weighted means (CWMs) along the dc-CA axes (composite gradients) and
 the matrix with linear combinations of environmental variables.
 \item \code{correlation}: a matrix with inter-set correlations of the environmental variables with their CWMs.
 \item \code{b_se}: a matrix with (unstandardized) regression coefficients for environmental
 variables and their optimistic standard errors.
 \item \code{R2_env}: a vector with coefficient of determination (R2) of the CWMs on to the environmental variables.
 The square-root thereof has been called the species-environmental correlation in CCA.
 }

}
All scores in the \code{dccav} object
are in scaling \code{"sites"} (1): the scaling with \emph{Focus on Case distances} .
}
\description{
Double constrained correspondence analysis (dc-CA) for analyzing
(multi-)trait (multi-)environment ecological data using library \code{vegan} and native R code.
It has a \code{formula} interface which allows to assess, for example, the importance
of trait interactions in shaping ecological communities.
The function \code{dc_CA} has an option to divide the abundance data of a site by the site total,
giving equal site weights.
This division has the advantage that the multivariate analysis corresponds with an unweighted (multi-trait) community-level analysis, instead of being weighted.
}
\details{
Empty (all zero) rows and columns in \code{response} are removed from the \code{response} and the corresponding
rows from \code{dataEnv} and \code{dataTraits}. Subsequently, any columns with missing values
are removed from  \code{dataEnv} and \code{dataTraits}. It gives an error (object 'name_of_variable' not found),
if variables with missing entries are specified in \code{formulaEnv} and \code{formulaTraits}.

Computationally, dc-CA can be carried out by a single singular value decomposition (ter Braak et al. 2018),
but it is here computed in two steps.
In the first step, the transpose of the \code{response}
is regressed on to the traits (the column predictors) using \code{\link[vegan]{cca}}
with \code{formulaTraits}.
The column scores of this analysis (in scaling 1) are community weigthed means (CWM) of the
orthonormalized traits.
These are then regressed on the environmental (row) predictors using \code{\link{wrda}} with
with \code{formulaEnv} or using \code{\link[vegan]{rda}}, if site weights are equal.

A dc-CA can be carried out on, what statisticians call, the sufficient statistics of the method.
In this case, \code{response} should be a list with as first element
community weighted means (CWMs) with respect to the traits, and the trait data,
and, optionally, further elements, for functions related to \code{dc_CA}.
The minimum is a \code{list(CWM, weight = list(columns= species_weights))}
with CWM a matrix or data.frame (or with data.frame in its class vector), but then
formulaEnv, formulaTraits, dataEnv, dataTraits must be specified in the call to \code{dc_CA}.
The function \code{\link{fCWM_SNC}} shows how to set the \code{response} for this and helps
to create the \code{response} from abundance data in these non-standard applications of dc-CA.
Species and site weights, if not set in \code{response$weights} can be set by a variable
\code{weight} in the data frames \code{dataTraits} and \code{dataEnv}, respectively,
but formulas should the not be \code{~.}.

The statistics and scores in the example \code{dune_dcCA.r},
 have been checked against the results in Canoco 5.15 (ter Braak & Smilauer, 1918).

In the current implementation, \code{formulaEnv} and \code{formulaTraits} should
contain variable names as is, \emph{i.e.} transformations of variables in the formulas gives
an error ('undefined columns selected') when the \code{\link{scores}} function is applied.
}
\examples{
library(douconca)
data("dune_trait_env")
# rownames are carried forward in results
rownames(dune_trait_env$comm) <- dune_trait_env$comm$Sites
divide <- FALSE # divide by site.totals if TRUE
cat("\n\n\n******* divide.by.site.totals ==", divide, "*******\n\n\n")
mod <- dc_CA(formulaEnv = ~A1+Moist+Mag+Use+Manure,
                   formulaTraits = ~. ,
                   response = dune_trait_env$comm[, -1],  # must delete "Sites"
                   dataEnv =dune_trait_env$envir,
# delete "Species", "Species_abbr" from traits and
# use all remaining variables due to formulaTraits = ~. (the default)
                   dataTraits =dune_trait_env$traits[,-c(1,2)],
                   divide.by.site.totals = divide,
                   verbose = TRUE)


class(mod)
print(mod)

set.seed(123)

anova(mod, by= "axis")
# For more demo on testing, see demo dune_test.r

# douconca:: was added explicitly as sometimes the vegan scores function takes over command.
# with error messaga:  'arg' should be one of "sites", "species", "both"
mod_scores <- douconca::scores(mod, display = c("all"), scaling = "sym")

douconca::scores(mod, display = c("cor", "cn", "cor_traits"),
                 scaling = "sym", which_cor = list("SLA","in model"))
douconca::scores(mod, display = c("cor", "cn", "cor_traits"),
                 scaling = "sym", which_cor = list("in model","Manure"))


cat("\n\n**meaning of each set of scores**\n")
print(meaning <- lapply(mod_scores, function(x)attr(x,which = "meaning")))

cat("head of unconstrained site scores, with meaning\n")
print(head(mod_scores$sites))
cat(attr(mod_scores$sites, which = "meaning"),"\n\n")
mod_scores$regression
mod_scores$biplot
mod_scores$centroids
# species related scores
cat("head of unconstrained species scores, with meaning\n")
print(head(mod_scores$species))
cat(attr(mod_scores$species, which = "meaning"),"\n\n")
mod_scores$regression_traits
mod_scores$biplot_traits
mod_scores$centroids_traits


mod_scores_tidy <- douconca::scores(mod, display = "all", tidy = TRUE)
print("names of the tidy scores")
print(names(mod_scores_tidy))
cat("\nThe levels of the tidy scores\n")
print(levels(mod_scores_tidy$score))


cat("\nFor illustration: a dc-CA model with a trait covariate\n")
mod2 <- dc_CA(formulaEnv = ~A1+Moist+Mag+Use+Manure,
                   formulaTraits = ~ SLA+Height+ LDMC+ Lifespan +Condition(Seedmass) ,
                   response = dune_trait_env$comm[, -1],  # must delete "Sites"
                   dataEnv =dune_trait_env$envir,
                   dataTraits =dune_trait_env$traits,
                   divide.by.site.totals = divide,
                   verbose = TRUE)

cat("\nFor illustration: a dc-CA model with both environmental and trait covariates\n")
mod3 <- dc_CA(formulaEnv = ~A1+Moist+Use+Manure+Condition(Mag),
                    formulaTraits = ~ SLA+Height+LDMC+Lifespan +Condition(Seedmass) ,
                    response = dune_trait_env$comm[, -1],  # must delete "Sites"
                    dataEnv =dune_trait_env$envir,
                    dataTraits =dune_trait_env$traits,
                    divide.by.site.totals = divide,
                    verbose = TRUE)

cat("\nFor illustration: same model but using mod2 for speed, as the trait model and data did not change\n")
mod3B <- dc_CA(formulaEnv = ~A1+Moist+Use+Manure+Condition(Mag),
                    dataEnv =dune_trait_env$envir,
                    dc_CA_object = mod2,
                    verbose = TRUE)
cat("\ncheck on equality of mod3 (from data) and mod3B (from a dc_CA_object)\n",
    "the expected difference is in the component 'call'\n ")
print(all.equal(mod3,mod3B)) #  an expected difference, namely in component call

mod_scores <- douconca::scores(mod3, display = "all", scaling = "symmetric")
#print(str(mod_scores))


cat("head of unconstrained site scores\n")
print(head(mod_scores$sites))
attr(mod_scores$sites, which = "meaning")
mod_scores$regression

mod_scores$biplot
mod_scores[["biplot"]]
mod_scores$centroids
# #Warning: gives the mod_scores$centroids_traits,
#           as mod_scores[["centroids"]] is NULL.
# Use mod_scores[["centroids"]] instead:
mod_scores[["centroids"]]
# species related scores
cat("head of unconstrained species scores\n")
print(head(mod_scores$species))
attr(mod_scores$species, which = "meaning")
mod_scores$regression_traits
mod_scores$biplot_traits
mod_scores$centroids_traits



mod_scores <- douconca::scores(mod3, display = "all", tidy = TRUE)
names(mod_scores)
levels(mod_scores$score)
cat("\nthe tidy scores data frame \n")
str(mod_scores)




# All statistics and scores have been checked against the results
# in Canoco 5.15 (ter Braak & Smilauer, 1918) for all three types of scaling.

}
\references{
Kleyer, M., Dray, S., Bello, F., Lepš, J., Pakeman, R.J., Strauss, B., Thuiller,
W. & Lavorel, S. (2012) Assessing species and community functional responses to
environmental gradients: which multivariate methods?
Journal of Vegetation Science, 23, 805-821.
http://dx.doi.org/10.1111/j.1654-1103.2012.01402.x

ter Braak, CJF, Šmilauer P, and Dray S. 2018. Algorithms and biplots for
double constrained correspondence analysis.
Environmental and Ecological Statistics, 25(2), 171-197.
https://doi.org/10.1007/s10651-017-0395-x or
http://rdcu.be/ETPh

ter Braak C.J.F. and  P. Šmilauer  (2018). Canoco reference manual
and user's guide: software for ordination (version 5.1x).
Microcomputer Power, Ithaca, USA, 536 pp.

Oksanen, J., et al. (2022)
vegan: Community Ecology Package. R package version 2.6-4.
http://CRAN.R-project.org/package=vegan.
}
\seealso{
\code{\link{plot_dcCA}}, \code{\link{scores.dcca}}, \code{\link{print.dcca}} and \code{\link{anova.dcca}}
}
