
<!-- README.md is generated from README.Rmd. Please edit that file -->

# douconca

<!-- badges: start -->
<!-- badges: end -->

R library `douconca` analyzes multi-trait multi-environment ecological
data by double constrained correspondence analysis (ter Braak et al.
2018) using `vegan` and native R code. It has a `formula` interface for
the trait- (column-) and environment- (row-) models, which allows to
assess, for example, the importance of trait interactions in shaping
ecological communities. Throughout the two step algorithm of ter Braak
et al. (2018) is used. This algorithm combines and extends community-
(sample-) and species-level analyses, i.e. the usual community weighted
means (CWM)-based regression analysis and the species-level analysis of
species-niche centroids (SNC)-based regression analysis. The CWM
regressions are specified with an environmental formula and the SNC
regressions are specified with a trait formula. dcCA finds the
environmental and trait gradients that optimize these regressions. The
first step uses
<a href="https://vegandevs.github.io/vegan/reference/cca.html"
class="uri">cca</a> (Oksanen et al. 2022) to regress the transposed
abundance data on to the traits and (weighted) redundancy analysis to
regress the community-weighted means (CWMs) of the orthonormalized
traits, obtained from the first step, on to the environmental
predictors. The sample total of the abundance data are used as weights.
The redundancy analysis is carried out using
<a href="https://vegandevs.github.io/vegan/reference/rda.html"
class="uri">rda</a> if sites have equal weights (after division of the
rows by their total) or, in the general weighted case, using `wrda`.
Division by the sample total has the advantage that the multivariate
analysis corresponds with an unweighted (multi-trait) community-level
analysis, instead of being weighted, which may give a puzzling
difference between common univariate and this multivariate analysis.

Reference: ter Braak, CJF, Šmilauer P, and Dray S. 2018. Algorithms and
biplots for double constrained correspondence analysis. Environmental
and Ecological Statistics, 25(2), 171-197.
<https://doi.org/10.1007/s10651-017-0395-x>

## Installation

You can install the development version of `douconca` like so:

``` r
install.packages("remotes")
remotes::install_github("CajoterBraak/douconca")
```
