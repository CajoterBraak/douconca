data("dune_trait_env")

# rownames are carried forward in results
rownames(dune_trait_env$comm) <- dune_trait_env$comm$Sites

mod <- dc_CA(formulaEnv = ~ A1 + Moist + Mag + Use + Manure,
             formulaTraits = ~ .,
             response = dune_trait_env$comm[, -1],  # must delete "Sites"
             dataEnv = dune_trait_env$envir,
             # delete "Species", "Species_abbr" from traits and
             # use all remaining variables due to formulaTraits = ~. (the default)
             dataTraits = dune_trait_env$traits[, -c(1,2)],
             verbose = TRUE)
anova(mod)

a_species <- anova_species(mod)
a_species
# anova_species can be used for manual forward selection of
# trait variables, as done for environmental variables in the demo 
# dune_FS_dcCA.r, based on the first eigenvalue and its significance
# and adding the first axis of the constrained species scores from mod to 
# the Condition of a new mod.
(eig1_and_pval <- c(eig = a_species$eig[1], pval = a_species$table$`Pr(>F)`[1]))
a_species$eig
anova_sites(mod)