library(douconca)
data("dune_trait_env")
# rownames are carried forward in results
rownames(dune_trait_env$comm) <- dune_trait_env$comm$Sites
divide <- FALSE # divide by site.totals if TRUE
cat("\n\n\n******* divide.by.site.totals ==", divide, "*******\n\n\n")
mod <- dc_CA(formulaEnv = ~A1+Moist+Mag+Use+Manure,
                   formulaTraits = ~. ,
                   response = dune_trait_env$comm[, -1],  # must delete "Sites"
                   dataEnv =dune_trait_env$envir,
# delete "Species", "Species_abbr" from traits and
# use all remaining variables due to formulaTraits = ~. (the default)
                   dataTraits =dune_trait_env$traits[,-c(1,2)],
                   divide.by.site.totals = divide,
                   verbose = TRUE)


class(mod)
print(mod)

set.seed(123)

anova(mod, by= "axis")
# For more demo on testing, see demo dune_test.r

# douconca:: was added explicitly as sometimes the vegan scores function takes over command.
# with error messaga:  'arg' should be one of "sites", "species", "both"
mod_scores <- douconca::scores(mod, display = c("all"), scaling = "sym")

douconca::scores(mod, display = c("cor", "cn", "cor_traits"),
                 scaling = "sym", which_cor = list("SLA","in model"))
douconca::scores(mod, display = c("cor", "cn", "cor_traits"),
                 scaling = "sym", which_cor = list("in model","Manure"))


cat("\n\n**meaning of each set of scores**\n")
print(meaning <- lapply(mod_scores, function(x)attr(x,which = "meaning")))

cat("head of unconstrained site scores, with meaning\n")
print(head(mod_scores$sites))
cat(attr(mod_scores$sites, which = "meaning"),"\n\n")
mod_scores$regression
mod_scores$biplot
mod_scores$centroids
# species related scores
cat("head of unconstrained species scores, with meaning\n")
print(head(mod_scores$species))
cat(attr(mod_scores$species, which = "meaning"),"\n\n")
mod_scores$regression_traits
mod_scores$biplot_traits
mod_scores$centroids_traits


mod_scores_tidy <- douconca::scores(mod, display = "all", tidy = TRUE)
print("names of the tidy scores")
print(names(mod_scores_tidy))
cat("\nThe levels of the tidy scores\n")
print(levels(mod_scores_tidy$score))


cat("\nFor illustration: a dc-CA model with a trait covariate\n")
mod2 <- dc_CA(formulaEnv = ~A1+Moist+Mag+Use+Manure,
                   formulaTraits = ~ SLA+Height+ LDMC+ Lifespan +Condition(Seedmass) ,
                   response = dune_trait_env$comm[, -1],  # must delete "Sites"
                   dataEnv =dune_trait_env$envir,
                   dataTraits =dune_trait_env$traits,
                   divide.by.site.totals = divide,
                   verbose = TRUE)

cat("\nFor illustration: a dc-CA model with both environmental and trait covariates\n")
mod3 <- dc_CA(formulaEnv = ~A1+Moist+Use+Manure+Condition(Mag),
                    formulaTraits = ~ SLA+Height+LDMC+Lifespan +Condition(Seedmass) ,
                    response = dune_trait_env$comm[, -1],  # must delete "Sites"
                    dataEnv =dune_trait_env$envir,
                    dataTraits =dune_trait_env$traits,
                    divide.by.site.totals = divide,
                    verbose = TRUE)

cat("\nFor illustration: same model but using mod2 for speed, as the trait model and data did not change\n")
mod3B <- dc_CA(formulaEnv = ~A1+Moist+Use+Manure+Condition(Mag),
                    dataEnv =dune_trait_env$envir,
                    dc_CA_object = mod2,
                    verbose = TRUE)
cat("\ncheck on equality of mod3 (from data) and mod3B (from a dc_CA_object)\n",
    "the expected difference is in the component 'call'\n ")
print(all.equal(mod3,mod3B)) #  an expected difference, namely in component call

mod_scores <- douconca::scores(mod3, display = "all", scaling = "symmetric")
#print(str(mod_scores))


cat("head of unconstrained site scores\n")
print(head(mod_scores$sites))
attr(mod_scores$sites, which = "meaning")
mod_scores$regression

mod_scores$biplot
mod_scores[["biplot"]]
mod_scores$centroids
# #Warning: gives the mod_scores$centroids_traits,
#           as mod_scores[["centroids"]] is NULL.
# Use mod_scores[["centroids"]] instead:
mod_scores[["centroids"]]
# species related scores
cat("head of unconstrained species scores\n")
print(head(mod_scores$species))
attr(mod_scores$species, which = "meaning")
mod_scores$regression_traits
mod_scores$biplot_traits
mod_scores$centroids_traits



mod_scores <- douconca::scores(mod3, display = "all", tidy = TRUE)
names(mod_scores)
levels(mod_scores$score)
cat("\nthe tidy scores data frame \n")
str(mod_scores)




# All statistics and scores have been checked against the results
# in Canoco 5.15 (ter Braak & Smilauer, 1918) for all three types of scaling.

