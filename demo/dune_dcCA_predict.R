data("dune_trait_env")

# rownames are carried forward in results
rownames(dune_trait_env$comm) <- dune_trait_env$comm$Sites

mod <- dc_CA(formulaEnv = ~ A1 + Moist + Mag + Use + Condition(Manure),
             formulaTraits = ~ SLA + Height + LDMC + Condition(Seedmass) + Lifespan,
             response = dune_trait_env$comm[, -1],  # must delete "Sites"
             dataEnv = dune_trait_env$envir,
             dataTraits = dune_trait_env$traits)

# regression coefficients
predict(mod, type = "reg_env")
predict(mod, type = "reg_traits")

# fit the mean traits at each site (20x6),
# that is CWM at each site
pred.traits <- predict(mod, type = "traits")
head(pred.traits)

# fit the mean environment for each species (28x8)
# that is SNC of each species
pred.env <- predict(mod, type = "env")
head(pred.env)

pred.resp <- predict(mod, type = "response")
# pred has negative values and dc_CA cannot have negatives in the response
# so, modify pred.resp,
#whichgives about similar eigenvalues as the original data
pred.resp[pred.resp < 0] <- 0
mod3 <- dc_CA(formulaEnv = mod$formulaEnv,
              formulaTraits = mod$formulaTraits,
              response = pred.resp, 
              dataEnv = dune_trait_env$envir,
              dataTraits = dune_trait_env$traits)
mod3$eigenvalues / mod$eigenvalues

