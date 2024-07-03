data("dune_trait_env")

# rownames are carried forward in results
rownames(dune_trait_env$comm) <- dune_trait_env$comm$Sites
response <- dune_trait_env$comm[, -1]  # must delete "Sites"

w <- rep(1, 20) 
w[1:10] <- 8 
w[17:20] <- 0.5

object <- wrda(formula = ~ A1 + Moist + Mag + Use + Condition(Manure),
               response = response, 
               data = dune_trait_env$envir, 
               weights = w)
object # Proportions equal to those Canoco 5.15

mod_scores <- scores(object, display = "all")
scores(object, which_cor = c("A1", "Manure"), display = "cor")
anova(object)

# The default is equal weights, which allows checking against vegan
object1 <- wrda(formula = ~ A1 + Moist + Mag + Use + Condition(Manure),
                response = response, 
                data = dune_trait_env$envir)
object1

# compare with vegan::rda
object2 <- 
  vegan::rda(formula = response ~ A1 + Moist + Mag + Use + Condition(Manure), 
             data = dune_trait_env$envir)
object2

object1$CCA$eig / object2$CCA$eig 
object1$tot.chi / object2$tot.chi
object1$CCA$tot.chi / object2$CCA$tot.chi
range(object1$CCA$u / object2$CCA$u)

myconst <- sqrt(object1$Nobs * object2$tot.chi)
veganScores <- 
  scores(object2, 
         choices = seq_len(ncol(object1$site_axes$site_scores$lc_env_scores)),
         display = "lc", scaling = "sites", const = myconst)
range(object1$site_axes$site_scores$lc_env_scores / veganScores) 