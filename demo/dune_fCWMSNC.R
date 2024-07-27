data("dune_trait_env")

# rownames are carried forward in results
rownames(dune_trait_env$comm) <- dune_trait_env$comm$Sites

CWMSNC <- fCWM_SNC(formulaEnv = ~ A1 + Moist + Manure + Use + Condition(Mag),
                   formulaTraits = ~ SLA + Height + LDMC + Condition(Seedmass) + Lifespan,
                   response = dune_trait_env$comm[, -1],  # must delete "Sites"
                   dataEnv = dune_trait_env$envir,
                   dataTraits = dune_trait_env$traits)
names(CWMSNC)
#CWMSNC$SNC <- NULL # would give correct dc-CA but no species-level t-values or test
mod <- dc_CA(response = CWMSNC) # formulas and data are in CWMSNC!
# note that output also gives the environment-constrained inertia,
# which is a bonus compare to the usual way to carry out a dcCA.
anova(mod)

