# douconca

<!-- badges: start -->
<!-- badges: end -->

 douconca analyzes multi-trait multi-environment ecological data by
double constrained correspondence analysis (ter Braak et al. 2018) 
using \code{vegan} and native R code. It has a \code{formula} interface
for the trait- (column-) and environment- (row-) models,
which allows to assess, for example, the importance
of trait interactions in shaping ecological communities.
Throughout the two step algorithm of ter Braak et al. (2018) is used. This algorithm
combines and extends community- (sample-) and species-level analyses, i.e.
the usual community weighted means (CWM)-based regression analysis and the
species-level analysis of species-niche centroids (SNC)-based regression analysis.
The CWM regressions are specified with an environmental formula 
and the SNC regressions are specified with a trait formula. dcCA finds 
the environmental and trait gradients that optimize these regressions.
The first step of the algorithm uses \code{\link[vegan]{cca}} 
to regress the (transposed) abundance data on to the traits 
and the second step uses \code{\link[vegan]{rda}} to regress the CWMs of the orthonormalized traits on to the environmental predictors.
The abundance data are divided by the sample total
(i.e. 'closed') in the vegan-based version. This
has the advantage that this multivariate analysis corresponds with an unweighted (multi-trait)
community-level analysis, instead of being weighted, which may give a puzzling differences
between common univarite and this multivariate analysis.
Reference: ter Braak, CJF, Šmilauer P, and Dray S. 2018. Algorithms and biplots fordouble constrained correspondence analysis. Environmental and Ecological Statistics, 25(2), 171-197. https://doi.org/10.1007/s10651-017-0395-x

## Installation

You can install the released version of dcCA from github by
invoking the R-code within an R-console:

``` r
install.packages("remotes")
remotes::install_github("CajoterBraak/douconca")
```

