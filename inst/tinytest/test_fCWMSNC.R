data("dune_trait_env")

# rownames are carried forward in results
rownames(dune_trait_env$comm) <- dune_trait_env$comm$Sites
divide = TRUE

dcca_mod <- dc_CA(formulaEnv = ~ A1 + Moist + Manure + Use + Condition(Mag),
                  formulaTraits = ~ SLA + Height + LDMC + Condition(Seedmass) + Lifespan,
                  response = dune_trait_env$comm[, -1],  # must delete "Sites"
                  dataEnv = dune_trait_env$envir,
                  dataTraits = dune_trait_env$traits,
                  divideBySiteTotals = divide,
                  verbose = FALSE)

CWMSNCa <- fCWM_SNC(formulaEnv = dcca_mod$formulaEnv,
                    formulaTraits = dcca_mod$formulaTraits,
                    response = dune_trait_env$comm[, -1],  # must delete "Sites"
                    dataEnv = dune_trait_env$envir,
                    dataTraits = dune_trait_env$traits,
                    divideBySiteTotals = divide,
                    verbose = FALSE)

expect_inherits(CWMSNCa, "list")
expect_equal(names(CWMSNCa), 
             c("CWM", "SNC", "formulaEnv", "formulaTraits", "env_explain", 
               "weights", "call", "data"))

CWMSNCb <- fCWM_SNC(formulaEnv = dcca_mod$formulaEnv,
                    formulaTraits = dcca_mod$formulaTraits,
                    response = dune_trait_env$comm[, -1],  # must delete "Sites"
                    dataEnv = dune_trait_env$envir,
                    dataTraits = dune_trait_env$traits,
                    divideBySiteTotals = divide,
                    minimal_output = FALSE,
                    verbose = FALSE)

expect_inherits(CWMSNCb, "list")
expect_equal(names(CWMSNCb), 
             c("CWM", "CWMs_orthonormal_traits", "SNC", "SNCs_orthonormal_env", 
               "Nobs", "traits_explain", "env_explain", "formulaEnv", 
               "formulaTraits", "trans2ortho", "T_ortho", "weights", "call", 
               "data"))

dcca_mod2 <- dc_CA(formulaEnv = dcca_mod$formulaEnv,
                   response = CWMSNCa,  
                   verbose = FALSE)

# Note: the axes of dcca_mod2 may have switched sign compared to dcca_mod
# formulaTrait is taken from CWMSNCa
CWMSNCc <- CWMSNCa
CWMSNCc$SNC <- NULL
dcca_mod3 <- dc_CA(formulaEnv = dcca_mod$formulaEnv,
                   response = CWMSNCc,  
                   verbose = FALSE)

expect_equal(dcca_mod$eigenvalues, dcca_mod2$eigenvalues)
expect_equal(dcca_mod$eigenvalues, dcca_mod3$eigenvalues)

expect_equal(abs(scores(dcca_mod2, display = "sites")),
             abs(scores(dcca_mod, display = "sites")))
expect_equal(abs(scores(dcca_mod2, display = "species")),
             abs(scores(dcca_mod, display = "species")))

# example of no weights specified
CWMSNCd <- CWMSNCb
CWMSNCd$weights <- NULL
expect_warning(dcca_mod4 <- dc_CA(formulaEnv = dcca_mod$formulaEnv,
                                  response = CWMSNCd, 
                                  verbose = FALSE),
               "no weights supplied")

expect_equivalent(as.numeric(dcca_mod4$eigenvalues), 
                  c(0.0549873463880468, 0.0241323582517041, 
                    0.0128375843092032, 0.0016183457949121))

# example of weights specified via dataTraits and dataEnv
CWMSNCe <- CWMSNCd
CWMSNCe$dataTraits <- dune_trait_env$traits
CWMSNCe$dataTraits$weight <- CWMSNCa$weights$columns
dune_trait_env$envir$weight <- CWMSNCa$weights$rows

expect_warning(dcca_mod5 <- dc_CA(formulaEnv = dcca_mod$formulaEnv,
                                  dataEnv = dune_trait_env$envir,
                                  response = CWMSNCe,  
                                  verbose = FALSE),
               "species weights taken from")

expect_equal(dcca_mod5$eigenvalues, dcca_mod$eigenvalues)

# A minimal specification with a non-trivial trait model, giving 3 warnings
dune_trait_env$envir$weight <- NULL
CWMSNCf <- list(CWM = as.data.frame(CWMSNCa$CWM),
                weights = list(columns = 100 * dcca_mod$weights$columns),
                dataTraits = dune_trait_env$traits,
                formulaTraits = ~ SLA + Height + LDMC + Condition(Seedmass) + Lifespan)

# Without trait covariates and only traits in dataTraits, formulaTraits can be deleted from the list.
# For the trait model, different environmental predictors can be used
expect_warning(dcca_mod6 <- dc_CA(response = CWMSNCf,
                                  dataEnv = dune_trait_env$envir,
                                  formulaEnv = ~ Moist,
                                  verbose = FALSE),
               "no site weights supplied")

