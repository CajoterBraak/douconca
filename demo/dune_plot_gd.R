data("dune_trait_env")

# rownames are carried forward in results
rownames(dune_trait_env$comm) <- dune_trait_env$comm$Sites

mod <- dc_CA(formulaEnv = ~A1 + Moist + Mag + Use + Manure,
             formulaTraits = ~ SLA + Height + LDMC + Seedmass + Lifespan,
             response = dune_trait_env$comm[, -1],  # must delete "Sites"
             dataEnv = dune_trait_env$envir,
             dataTraits = dune_trait_env$traits,
             verbose = FALSE)

env_scores <- scores(mod, display = "tval")

env_scores <- data.frame(env_scores)
env_scores$group <- c("quantitative", "category")[c(1,1,2,2,2,1,1)]

plot_species_scores_bk(
  species_scores = env_scores,
  ylab = "optimistic t-values",  threshold = 0,  y_lab_interval = 1,
  scoresname = "dcCA1", speciesgroup = "group", verbose = FALSE
)
