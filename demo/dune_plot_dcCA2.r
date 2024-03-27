
data("dune_trait_env")
# rownames are carried forward in results
rownames(dune_trait_env$comm) <- dune_trait_env$comm$Sites
# must delete "Sites" from response matrix or data frame
Y <- dune_trait_env$comm[,-1] # must delete "Sites"

out <- dc_CA_vegan(formulaEnv = ~A1+Moist+Use+Manure+Condition(Mag),
                   formulaTraits = ~ SLA + Height + LDMC + Seedmass +Lifespan,
                   response = Y,
                   dataEnv =dune_trait_env$envir,
                   dataTraits =dune_trait_env$traits,
                   verbose = FALSE)

#print(out)

plot_dcCA_CWM_SNC(out, facet = FALSE)
CWM_SNC_env_trait_scores <- getPlotdata(out)
myplot<-plot_dcCA(out, verbose = FALSE)
names(myplot)
plot_dcCA(out, gradient_description = c("corre", "tval"), verbose = FALSE )

newnames_without_covariates <-list(traits= c("SLA", "Height", "LDMC", "Seedmass", "annual", "perennial"),
                env= c("A1 horizon", "Moisture", "Type of use",  "Manure"))

# modifying the plot
# Assign the plots to symbols and use grid.arrange to produce the plot you like, for example:
gg<- plot_dcCA(out, gradient_description = c("corre", "weights"), verbose = FALSE, newnames = newnames_without_covariates)
names(gg)
pl <- gg$separateplots
names(pl)
layout<- rbind(c(1,2,4), c(1,3,4))
gg_object <- gridExtra::grid.arrange(pl$CWM_SNC + ggplot2::xlab("env. and trait gradient"),
                                     pl$traits + ggplot2::ggtitle("cor"),
                                     pl$env+  ggplot2::ggtitle("w"),
                                     pl$species + ggplot2::ylab("trait complex"),
                                     layout_matrix = layout, widths = c(3,1,2),
                                     top = "",
                                     left ="", right =  "")




