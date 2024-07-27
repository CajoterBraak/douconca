data("dune_trait_env")

# rownames are carried forward in results
rownames(dune_trait_env$comm) <- dune_trait_env$comm$Sites

# must delete "Sites" from response matrix or data frame
Y <- dune_trait_env$comm[, -1] # must delete "Sites"

out <- dc_CA(formulaEnv = ~ A1 + Moist + Use + Manure + Mag,
                   formulaTraits = ~ SLA + Height + LDMC + Seedmass + Lifespan,
                   response = Y,
                   dataEnv = dune_trait_env$envir,
                   dataTraits = dune_trait_env$traits,
                   verbose = FALSE)
dat <- getPlotdata(out)
names(dat)
names(dat$CWM_SNC)
levels(dat$CWM_SNC$groups)

plot(out, verbose = FALSE)
