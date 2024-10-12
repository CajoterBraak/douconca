# douconca 1.2.2

* New coef.dcca() and fitted.dcca() functions with predict.dcca() adapted.
The function coef() can give fourth-corner correlations and regression 
coefficients. 
* Patch release with extended test files and associated small corrections.

# douconca 1.2.1

* Patch release addressing check errors on several CRAN build machines.

# douconca 1.2.0

* Intial CRAN release

# douconca 1.1.6

* An issue with collinear predictors in v1.1.5 has been resolved.

# douconca 1.1.5

* The package can now do general dc-CA, instead of the vegan-based version with 
equal site weights only. For users of the previous version, the function 
dc_CA_vegan has been replaced by the more general function dc_CA. 
The default gives the same analysis. By specifying
the argument `divideBySiteTotals = FALSE`, obtain the original dc-CA analysis
with unequal site weights.
* The `plot_dcCA` function is now a method: `plot.`
* General dc-CA required weighted redundancy analysis. For this, a new function
`wrda` has been added, with methods for print, scores and anova.
* A `predict` function has been added.
* A dc-CA can be computed from community-weighted means (CWMs) with
trait and environment data with species and site weights. See the new function 
`fCWM_SNC`. This is of interest, for example, to make a dc-CA analysis 
reproducible when the abundance data cannot be made public, and
it may also allow to perform dcCA with intra-species trait variation. 
The user needs to be able to compute meaningful CWMs in this case and supply 
trait data that reflect the (species-weighted) inter-trait covariance.
* Several functions are updated. In particular, there are corrections to
the anova function.

# douconca 1.1.2

* The `scores.dccav` function is corrected concerning intra-set correlations for
traits and environmental variables.
* The plotting functions are updated to avoid ggplot2 warnings on color and 
size.
* The fitted straight lines in the plots use the implicit weights 
(they did already, but the help said they did not).

