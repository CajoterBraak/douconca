data("dune_trait_env")

# rownames are carried forward in results
rownames(dune_trait_env$comm) <- dune_trait_env$comm$Sites
response <- dune_trait_env$comm[, -1]  

w <- rep(1, 20)
w[1:10] <- 4 
w[17:20] <- 0.5

mod_wrda <- wrda(formula = ~ A1 + Moist + Manure + Use + Condition(Mag),
                 response = response, 
                 data = dune_trait_env$envir, 
                 weights = w)

expect_inherits(mod_wrda, "wrda")

expect_equivalent(mod_wrda$eig, anova(mod_wrda)$eig)

expect_equal_to_reference(mod_wrda, "mod_wrda")

expect_equal(anova(mod_wrda, by = "axis")$table$Variance,
             c(13.4759764625272, 3.51996709660156, 2.83281997217771, 
               1.90283508422467, 24.4005894936307))

scores_wrda <- scores(mod_wrda, display = "all")
expect_equal_to_reference(scores_wrda, "scores_wrda")

scores_sub <- scores(mod_wrda, which_cor = c("A1", "Manure"), display = "cor")
expect_equivalent(scores_wrda$correlation[c("A1", "Manure"), ], scores_sub)

# The default is equal weights, which allows checking against vegan
mod_wrda_ew <- wrda(formula = ~ A1 + Moist + Mag + Use + Condition(Manure),
                    response = response, 
                    data = dune_trait_env$envir)

# compare with vegan::rda
mod_vegan <- vegan::rda(formula = response ~ A1 + Moist + Mag + Use + Condition(Manure), 
                        data = dune_trait_env$envir)

expect_equivalent(mod_wrda_ew$CCA$eig, mod_vegan$CCA$eig)
expect_equal(mod_wrda_ew$tot.chi, mod_vegan$tot.chi)
expect_equal(mod_wrda_ew$CCA$tot.chi, mod_vegan$CCA$tot.chi)
expect_equivalent(mod_wrda_ew$CCA$u, mod_vegan$CCA$u)

const <- sqrt(mod_wrda_ew$Nobs * mod_vegan$tot.chi)
expect_equivalent(mod_wrda_ew$site_axes$site_scores$lc_env_scores,
                  scores(mod_vegan, 
                         choices = seq_len(ncol(mod_wrda_ew$site_axes$site_scores$lc_env_scores)),
                         display= "lc", scaling = "sites", const = const))

wrda_print <- print(mod_wrda)
expect_equal(names(wrda_print), 
             c("call", "method", "tot.chi", "formula", "site_axes", 
               "species_axes", "Nobs", "eigenvalues", "weights", "data", "eY", 
               "pCCA", "CCA", "CA", "c_env_normed"))
