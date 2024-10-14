## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(douconca)

## -----------------------------------------------------------------------------
library(douconca)
data("dune_trait_env")
names(dune_trait_env)
dim(dune_trait_env$comm[, -1]) ## without the variable "Sites"
dim(dune_trait_env$traits)
dim(dune_trait_env$envir)
names(dune_trait_env$traits)
names(dune_trait_env$envir)

## -----------------------------------------------------------------------------
Y <- dune_trait_env$comm[, -1] # must delete "Sites"
mod <- dc_CA(formulaEnv = ~ A1 + Moist + Use + Manure + Mag,
             formulaTraits = ~ SLA + Height + LDMC + Seedmass + Lifespan,
             response = Y,
             dataEnv = dune_trait_env$envir,
             dataTraits = dune_trait_env$traits)

## -----------------------------------------------------------------------------
mod$eigenvalues

## ----fig.width=7--------------------------------------------------------------
plot(mod,gradient_description = "t")

## -----------------------------------------------------------------------------
set.seed(1)
anova(mod, by = "axis")

## -----------------------------------------------------------------------------
mod_e <-  dc_CA(formulaEnv = ~ A1 + Moist + Manure + Use + Mag,
                formulaTraits = ~ F + R + N + L,
                response = Y,
                dataEnv = dune_trait_env$envir,
                dataTraits = dune_trait_env$traits)

## -----------------------------------------------------------------------------
anova(mod_e, by = "axis")$max

## -----------------------------------------------------------------------------
round(sqrt(mod_e$eigenvalues[1]), 2)

## -----------------------------------------------------------------------------
mod_mGe <-  dc_CA(formulaEnv = ~ A1 + Moist + Manure + Use + Mag,
                 formulaTraits = 
                   ~ SLA + Height + LDMC + Seedmass + Lifespan + Condition(F+R+N+L),
                 response = Y,
                 dataEnv = dune_trait_env$envir,
                 dataTraits = dune_trait_env$traits, verbose = FALSE)
anova(mod_mGe, by= "axis")$max

## -----------------------------------------------------------------------------
mod_LDMC <- dc_CA(formulaEnv = ~ A1 + Moist + Manure + Use + Mag,
                   formulaTraits = ~ LDMC,
                   response = Y, 
                   dataEnv = dune_trait_env$envir,
                   dataTraits = dune_trait_env$trait, verbose = FALSE)
anova(mod_LDMC)

## -----------------------------------------------------------------------------
CWMSNC_LDMC <- fCWM_SNC(formulaEnv = ~ A1 + Moist + Manure + Use + Mag,
                      formulaTraits = ~ LDMC,
                      response = Y, 
                      dataEnv = dune_trait_env$envir,
                      dataTraits = dune_trait_env$trait, verbose = FALSE)

## -----------------------------------------------------------------------------
envCWM <- cbind(dune_trait_env$envir, CWMSNC_LDMC$CWM)
lmLDMC <- lm(LDMC ~ A1 + Moist + Manure + Use + Mag, data = envCWM)
anova(lmLDMC, lm(LDMC ~ 1, data = envCWM))

## -----------------------------------------------------------------------------
(regr_table <- scores(mod_LDMC, display = "reg"))
coefs_LDMC_dcCA <- regr_table[, "dcCA1"] / regr_table[, "SDS"]
range(coef(lmLDMC)[-1] / coefs_LDMC_dcCA)

## -----------------------------------------------------------------------------
cbind(summary(lmLDMC)$coefficients[-1, "t value", drop = FALSE],
scores(mod_LDMC, display = "tval"))

