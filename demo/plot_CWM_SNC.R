
data("dune_trait_env")
# rownames are carried forward in results
rownames(dune_trait_env$comm) <- dune_trait_env$comm$Sites
# must delete "Sites" from response matrix or data frame
Y <- dune_trait_env$comm[,-1] # must delete "Sites"

out <- dc_CA_vegan(formulaEnv = ~A1+Moist+Use+Manure+Condition(Mag),
                   formulaTraits = ~ SLA + Height + LDMC + Seedmass +Lifespan,
                   response = Y,
                   dataEnv =dune_trait_env$envir,
                   dataTraits =dune_trait_env$traits,
                   verbose = FALSE)
plot_dcCA_CWM_SNC(out, facet = FALSE)

newnames_without_covariates <-list(traits= c("SLA", "Height", "LDMC", "Seedmass", "annual", "perennial"),
                                   env= c("A1 horizon", "Moisture", "Type of use",  "Manure"))
plot_dcCA_CWM_SNC(out, newnames= newnames_without_covariates)

CWM_SNC_env_trait_scores <- getPlotdata(out)
