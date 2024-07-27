data("dune_trait_env")

# rownames are carried forward in results
rownames(dune_trait_env$comm) <- dune_trait_env$comm$Sites
divide <- TRUE # divide by site.totals if TRUE

mod_dcca <- dc_CA(formulaEnv = ~A1 + Moist + Mag + Use + Manure,
                  formulaTraits = ~ SLA + Height + LDMC + Seedmass + Lifespan,
                  response = dune_trait_env$comm[, -1],  # must delete "Sites"
                  dataEnv = dune_trait_env$envir,
                  dataTraits = dune_trait_env$traits,
                  divideBySiteTotals = divide,
                  verbose = FALSE)

predTraits <- predict(mod_dcca, type = "traits")
expect_equal_to_reference(predTraits, "predict_traits_dcca")

predEnv <- predict(mod_dcca, type = "env")
expect_equal_to_reference(predEnv, "predict_env_dcca")

predRegTraits <- predict(mod_dcca, type = "reg_traits")
expect_equal_to_reference(predRegTraits, "predict_regTraits_dcca")

predRegEnv <- predict(mod_dcca, type = "reg_env")
expect_equal_to_reference(predRegEnv, "predict_regEnv_dcca")

predResponse <- predict(mod_dcca, type = "response")
expect_equal_to_reference(predResponse, "predict_response_dcca")

CWMSNC <- list(CWM = predTraits, SNC = predEnv, weights = mod_dcca$weights, 
               data = mod_dcca$data, formulaEnv = mod_dcca$formulaEnv, 
               formulaTraits = mod_dcca$formulaTraits)
CWMSNC$data$Y <- NULL

expect_warning(mod_dcca2 <- dc_CA(response = CWMSNC,
                                  divideBySiteTotals = divide,
                                  verbose = FALSE),
               "overfitted model")

expect_equal(mod_dcca2$eigenvalues, mod_dcca$eigenvalues)
expect_null(mod_dcca2$site_axes$R2_env) # perfect fit
expect_null(mod_dcca2$species_axes$R2_traits) # perfect fit

# some variable may not be present in the newdata
# delete variable "Moist"
expect_silent(predict(mod_dcca, type = "traits", 
                      newdata = mod_dcca$data$dataEnv[, -3]))


predTraitsND <- predict(mod_dcca, type =  "traits", 
                        newdata = mod_dcca$data$dataTraits)
expect_equal(predTraitsND[1, rownames(mod_dcca$c_traits_normed)],
             mod_dcca$c_traits_normed[, "Avg"])

predEnvND <- predict(mod_dcca, type = "env", newdata  = mod_dcca$data$dataEnv)
expect_equivalent(predEnvND[1, ], c(4.85, 2.9, 0.3, 0.15, 0.25, 0.3, 1.9, 1.75))


## Leave out Moisture
mod_dcca2 <- dc_CA(formulaEnv = ~A1 + Mag + Use + Manure,
                   formulaTraits = ~ Height + LDMC +Seedmass +Lifespan,
                   response = dune_trait_env$comm[, -1],  # must delete "Sites"
                   dataEnv = dune_trait_env$envir,
                   dataTraits = dune_trait_env$traits,
                   divideBySiteTotals = divide,
                   verbose = FALSE)


predTraits2 <- predict(mod_dcca, type =  "traits", 
                       newdata  = mod_dcca$data$dataEnv[,-3]) # wo Moist
predEnv2 <- predict(mod_dcca, type =  "env", 
                    newdata  = mod_dcca$data$dataTraits[,-3]) # wo SLA
expect_equal(colnames(predTraits), 
             c("SLA", "Height", "LDMC", "Seedmass", 
               "Lifespanannual", "Lifespanperennial"))
expect_equal(colnames(predEnv2),
             c("A1", "Moist", "MagSF", "MagBF", "MagHF", "MagNM", 
               "Use", "Manure"))

predTraits2b <- predict(mod_dcca, type =  "traits", 
                        newdata  = mod_dcca$data$dataEnv[, 4, drop = FALSE]) # Mag only
expect_equal(nrow(unique(predTraits2b)), 4) # the four classes of Mag

predEnv2b <- predict(mod_dcca, type = "env", 
                     newdata = mod_dcca$data$dataTraits[, 7, drop = FALSE]) # LifeSpan only
expect_equal(nrow(unique(predEnv2b)), 2)  # the two classes of LifeSpan

CWMSNC3 <- list(CWM = predTraits2[, -1], SNC = predEnv2[, -2], 
                weights = mod_dcca$weights, data = mod_dcca$data,
                formulaEnv = mod_dcca2$formulaEnv, 
                formulaTraits = mod_dcca2$formulaTraits)
CWMSNC3$data$Y <- NULL

expect_warning(mod_dcca3 <- dc_CA(response = CWMSNC3,
                                  divideBySiteTotals = divide,
                                  verbose = FALSE),
               "overfitted model")

# the contribution of Moisture
expect_equal(mod_dcca$inertia["constraintsTE",1] - 
               mod_dcca3$inertia["constraintsTE", 1],
             0.295896752725216)

expect_equal(mod_dcca2$inertia["constraintsTE", 1] - 
               mod_dcca$inertia["constraintsTE", 1],
             -0.239836626665131)
