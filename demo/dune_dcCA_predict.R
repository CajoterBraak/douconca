data("dune_trait_env")

# rownames are carried forward in results
rownames(dune_trait_env$comm) <- dune_trait_env$comm$Sites

divide <- TRUE # divide by site.totals if TRUE

cat("\n\n\n******* divide.by.site.totals ==", divide, "*******\n\n\n")
mod <- dc_CA(formulaEnv = ~A1 + Moist + Mag + Use + Condition(Manure),
             formulaTraits = ~ SLA + Height + LDMC + Condition(Seedmass) + Lifespan,
             response = dune_trait_env$comm[, -1],  # must delete "Sites"
             dataEnv = dune_trait_env$envir,
             dataTraits = dune_trait_env$traits,
             divide.by.site.totals = divide,
             verbose = TRUE)

# regression coefficients
predict(mod, type = "reg_env")
predict(mod, type = "reg_traits")

# fit the mean traits at each site (20x6),
# that is CWM at each site, so analyse with CWMSNC version
pred.traits <- predict(mod, type = "traits")
head(pred.traits)

# fit the mean environment for each species (28x8)
# that is SNC of each species
pred.env <- predict(mod, type = "env")
head(pred.env)

CWMSNC <- list(CWM = pred.traits, SNC = pred.env, weights = mod$weights, 
               data = mod$data, formulaEnv = mod$formulaEnv, 
               formulaTraits = mod$formulaTraits)
CWMSNC$data$Y <- NULL

mod2 <- dc_CA(response = CWMSNC,
              divide.by.site.totals = divide,
              verbose = TRUE)

mod$eigenvalues / mod2$eigenvalues
mod2$site_axes$R2_env # perfect fit
mod2$species_axes$R2_traits # perfect fit

pred <- predict(mod, type = "response")
# pred has negative values and dc_CA cannot have negatives in the response
pred[pred < 0] <- 0

mod3 <- dc_CA(formulaEnv = mod$formulaEnv,
              formulaTraits = mod$formulaTraits,
              response = pred,  # must delete "Sites"
              dataEnv = dune_trait_env$envir,
              dataTraits = dune_trait_env$traits,
              divide.by.site.totals = divide,
              verbose = TRUE)

mod3$eigenvalues / mod$eigenvalues

