
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
class(mod)
print(mod)
set.seed(123)
# overall test
# community-level permutation test
anova(mod$RDAonEnv) # all options of anova.cca are available!
# by axis test
p_sites   <- anova(mod$RDAonEnv, by = "axis")
# a species-level permutation test required an dedicated new function
# anova_species(mod) # see dune_test.r

mod_scores <- scores(mod, display = c("all"), scaling = "sym")


(meaning <- lapply(mod_scores, function(x)attr(x,which = "meaning")))

head(mod_scores$sites)
attr(mod_scores$sites, which = "meaning")
mod_scores$regression
mod_scores$biplot
mod_scores$centroids
# species related scores
head(mod_scores$species)
attr(mod_scores$species, which = "meaning")
mod_scores$regression_traits
mod_scores$biplot_traits
mod_scores$centroids_traits


mod_scores_tidy <- scores(mod, display = "all", tidy = TRUE)
names(mod_scores_tidy)
levels(mod_scores_tidy$score)


# for illustration: a dc-CA model with a trait covariate
mod2 <- dc_CA_vegan(formulaEnv = ~A1+Moist+Mag+Use+Manure,
                   formulaTraits = ~ SLA+Height+ LDMC+ Lifespan +Condition(Seedmass) ,
                   response = dune_trait_env$comm[, -1],  # must delete "Sites"
                   dataEnv =dune_trait_env$envir,
                   dataTraits =dune_trait_env$traits,
                   verbose = TRUE)

# for illustration: a dc-CA model with both environmental and trait covariates
mod3 <- dc_CA_vegan(formulaEnv = ~A1+Moist+Use+Manure+Condition(Mag),
                    formulaTraits = ~ SLA+Height+LDMC+Lifespan +Condition(Seedmass) ,
                    response = dune_trait_env$comm[, -1],  # must delete "Sites"
                    dataEnv =dune_trait_env$envir,
                    dataTraits =dune_trait_env$traits,
                    verbose = TRUE)

# for illustration: same model but using mod2 for speed, as the trait model and data did not change
mod3B <- dc_CA_vegan(formulaEnv = ~A1+Moist+Use+Manure+Condition(Mag),
                    dataEnv =dune_trait_env$envir,
                    dc_CA_vegan_object = mod2,
                    verbose = TRUE)
all.equal(mod3,mod3B) # TRUE

mod_scores <- scores(mod3, display = "all", scaling = "symmetric")
str(mod_scores)


head(mod_scores$sites)
attr(mod_scores$sites, which = "meaning")
mod_scores$regression

mod_scores[["biplot"]]
#  mod_scores$centroids
# #gives the mod_scores$centroids_traits
# as mod_scores[["centroids"]] is NULL
mod_scores[["centroids"]]
# species related scores
head(mod_scores$species)
attr(mod_scores$species, which = "meaning")
mod_scores$regression_traits
mod_scores$biplot_traits
mod_scores$centroids_traits



mod_scores <- scores(mod3, display = "all", tidy = TRUE)
names(mod_scores)
levels(mod_scores$score)
str(mod_scores)




# All statistics and scores have been checked against the results
# in Canoco 5.15 (ter Braak & Smilauer, 1918) for all three types of scaling.

