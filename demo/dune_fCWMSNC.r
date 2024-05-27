library(douconca)
data("dune_trait_env")
# rownames are carried forward in results
rownames(dune_trait_env$comm) <- dune_trait_env$comm$Sites
divide = TRUE
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
              response = CWMSNCa,  # included dataEnv, dataTraits
              verbose = TRUE)
#Note: the axes of mod2 may have switched sign compared to mod1
# formulaTrait is taken from CWMSNCa
CWMSNCc <- CWMSNCa
CWMSNCc$SNC <- NULL
mod3 <- dc_CA(formulaEnv =mod1$formulaEnv,
              response = CWMSNCc,  # included dataEnv, dataTraits
              verbose = TRUE)
cbind(eig_mod1= mod1$eigenvalues, eig_mod2 = mod2$eigenvalues, eig_mod3= mod3$eigenvalues)
mod2$eigenvalues/ mod1$eigenvalues
type_of_score <- "sites"
douconca::scores(mod2, display= type_of_score)/douconca::scores(mod1, display= type_of_score)
type_of_score <- "species"
douconca::scores(mod2, display= type_of_score)/douconca::scores(mod1, display= type_of_score)
#Note: the axes of mod2 have switched sign compared to mod1


# example of no weights specified;
CWMSNCd <-CWMSNCb
CWMSNCd$weights <- NULL
mod4 <- dc_CA(formulaEnv =mod1$formulaEnv,
              response = CWMSNCd,  # included dataEnv, dataTraits
              verbose = FALSE)
# the analysis is not a real dc-CA:
# the eigenvalues of CWM and SNC analyses differ, see above warning
# the eigenvalues of the CWM analyses are not those of mod1
mod4$eigenvalues/mod1$eigenvalues


CWMSNCe <-CWMSNCd
CWMSNCe$weights <- NULL
# example of weights specified via dataTraits and dataEnv
CWMSNCe$dataTraits <- dune_trait_env$traits
CWMSNCe$dataTraits$weight <- CWMSNCa$weights$columns
dune_trait_env$envir$weight <- CWMSNCa$weights$rows
mod5 <- dc_CA(formulaEnv =mod1$formulaEnv,
              dataEnv = dune_trait_env$envir,
              response = CWMSNCe,  # included dataEnv, dataTraits
              verbose = TRUE)
mod5$eigenvalues/mod1$eigenvalues # mod5 is a proper dc-CA

# a minimal specification giving 4 warnings
dune_trait_env$envir$weight <- NULL
CWMSNCf <- list(CWM = as.data.frame(CWMSNC0$CWM),
                weights= list(columns = 100*mod1$weights$columns),
                dataTraits =dune_trait_env$traits[,-c(1,2)]
)
# the trait model is implicit in the CWMs, different environmental predictors can be used
mod6 <- dc_CA(response =CWMSNCf,
                dataEnv =dune_trait_env$envir,
                formulaEnv = ~ Moist,
                verbose = TRUE)

