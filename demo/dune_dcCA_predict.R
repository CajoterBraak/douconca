data("dune_trait_env")

# rownames are carried forward in results
rownames(dune_trait_env$comm) <- dune_trait_env$comm$Sites

mod <- dc_CA(formulaEnv = ~ A1 + Moist + Mag + Use + Condition(Manure),
             formulaTraits = ~ SLA + Height + LDMC + Condition(Seedmass) + Lifespan,
             response = dune_trait_env$comm[, -1],  # must delete "Sites"
             dataEnv = dune_trait_env$envir,
             dataTraits = dune_trait_env$traits, verbose = FALSE)

# Ten 'new' sites with a subset of the variables in mod 
# X_lot will be ignored as it is not part of the model
newEnv <- dune_trait_env$envir[1:10,c("A1", "Mag", "Manure", "X_lot")]
newEnv[2,"A1"] <- 3.0; rownames(newEnv) <- paste0("NewSite", 1:10)
pred.traits <- predict(mod, type = "traitsFromEnv", newdata = newEnv)
head(pred.traits)

# Eight 'new' species with a subset of traits that are included in the model 
# Variable "L" will be ignored as it is not in the model 
newTraits <- dune_trait_env$traits[1:8,c("Species","SLA", "LDMC", "L")]
newTraits[3,"SLA"]<- 18; 
rownames(newTraits) <- paste("Species",LETTERS[1:8] )# or another meaningful name...
pred.env <- predict(mod, type = "envFromTraits", newdata = newTraits)
head(pred.env)

pred.resp <- predict(mod, type = "response", newdata= list(newTraits,newEnv),
        weights= list(species = rep(c(1,2),4) , sites = rep(1,10)) )
colSums(pred.resp) # about alternating 0.8 and 1.6 (reflecting the new species weights)
rowSums(pred.resp) # about equal rowsums


