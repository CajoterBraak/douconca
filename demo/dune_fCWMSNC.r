data("dune_trait_env")
# rownames are carried forward in results
rownames(dune_trait_env$comm) <- dune_trait_env$comm$Sites
divide = FALSE
mod1 <- dc_CA(formulaEnv = ~A1+Moist+Manure+Use+Condition(Mag),
              formulaTraits = ~ SLA + Height + LDMC + Condition(Seedmass) +Lifespan,
              response = dune_trait_env$comm[, -1],  # must delete "Sites"
              dataEnv =dune_trait_env$envir,
              dataTraits =dune_trait_env$traits,
              divide.by.site.totals = divide,
              verbose = TRUE)


CWMSNCa <- fCWM_SNC(formulaEnv = mod1$formulaEnv,
                          formulaTraits = mod1$formulaTraits,
                          response = dune_trait_env$comm[, -1],  # must delete "Sites"
                          dataEnv =dune_trait_env$envir,
                          dataTraits =dune_trait_env$traits,
                          divide.by.site.totals = divide,
                          verbose = FALSE)
names(CWMSNCa)

CWMSNCb <- fCWM_SNC(formulaEnv = mod1$formulaEnv,
                          formulaTraits = mod1$formulaTraits,
                          response = dune_trait_env$comm[, -1],  # must delete "Sites"
                          dataEnv =dune_trait_env$envir,
                          dataTraits =dune_trait_env$traits,
                          divide.by.site.totals = divide,
                          minimal_output = FALSE,
                          verbose = FALSE)
names(CWMSNCb)

mod2 <- dc_CA(formulaEnv =mod1$formulaEnv,
               formulaTraits = mod1$formulaTraits,
               response = CWMSNCa,  # included dataEnv, dataTraits
               verbose = TRUE)
#Note: the axes of mod2 may have switched sign compared to mod1
CWMSNCc <- CWMSNCa
CWMSNCc$SNC <- NULL
mod3 <- dc_CA(formulaEnv =mod1$formulaEnv,
              formulaTraits = mod1$formulaTraits,
              response = CWMSNCc,  # included dataEnv, dataTraits
              verbose = TRUE)
cbind(eig_mod1= mod1$eigenvalues, eig_mod2 = mod2$eigenvalues, eig_mod3= mod3$eigenvalues)
mod2$eigenvalues/ mod1$eigenvalues
type_of_score <- "sites"
douconca::scores(mod2, display= type_of_score)/douconca::scores(mod1, display= type_of_score)
type_of_score <- "species"
douconca::scores(mod2, display= type_of_score)/douconca::scores(mod1, display= type_of_score)
#Note: the axes of mod2 have switched sign compared to mod1
