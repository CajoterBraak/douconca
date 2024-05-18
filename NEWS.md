# douconca 1.1.5

* The package can now do general dc-CA, instead of the vegan-based version with equal
site weights only. For users of the previous version, the function 
dc_CA_vegan is now replaced by the more general function dc_CA. Specify
the argument divide.by.site.totals to obtain the same analysis.
* General dc-CA required weighted redundancy analysis. For this a new function
wrda has been added, with methods for print, scores and anova.
*A dc-CA can be computed from community-weighted means with
trait and environment data. See the new function fCWM_SNC. This is of interest, for example,
to make a dc-CA analysis reproducible when the abundance data cannot be made public.
* Several functions are updated. In particular, there are correction to
the anova function.

# douconca 1.1.2

* The scores.dccav function is corrected concerning intra-set correlations for traits and environmental variables.
* The plotting functions are updated to avoid ggplot2 warnings on color and size.
* The fitted straight lines in the plots use the implicit weights (they did already, but the help said they did not.)

