data("dune_trait_env")
# rownames are carried forward in results
rownames(dune_trait_env$comm) <- dune_trait_env$comm$Sites
mod <- dc_CA_vegan(formulaEnv = ~A1+Moist+Mag+Use+Manure,
                   formulaTraits = ~. ,
                   response = dune_trait_env$comm[, -1],  # must delete "Sites"
                   dataEnv =dune_trait_env$envir,
# delete "Species", "Species_abbr" from traits and
# use all remaining variables due to formulaTraits = ~. (the default)
                   dataTraits =dune_trait_env$traits[,-c(1,2)],
                   verbose = TRUE)
anova(mod)

a_species <- anova_species(mod)
a_species
# anova_species can be used for manual forward selection of
# trait variables, as done for environmental variables in the demo dune_FS_dcCA.r,
# based on the first eigenvalue and its significance
# and adding the first axis of the constrained species scores from mod to the Condition
# of a new mod.
(eig1_and_pval <- c(eig = a_species$eig[1], pval = a_species$table$`Pr(>F)`[1]))
a_species$eig
#
anova_sites(mod)

mod$eigenvalues / anova_species(mod)$eig[seq_along(mod$eigenvalues)]
mod$eigenvalues / anova_sites(mod)$eig[seq_along(mod$eigenvalues)]








# Benchmark the community-level test
#microbenchmark::microbenchmark(anova(mod$RDAonEnv, by = "axis"), anova_sites(mod, by ="axis"),times = 1)
#microbenchmark::microbenchmark(anova(mod$RDAonEnv), anova_sites(mod), times = 2)
