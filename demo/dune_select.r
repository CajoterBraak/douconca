
data("dune_trait_env")

rownames(dune_trait_env$comm) <- dune_trait_env$comm$Sites
out1 <- dc_CA_vegan(formulaEnv = ~A1+Moist+Mag+Use+Manure,
                   formulaTraits = ~. ,
                   response = dune_trait_env$comm[, -1],  # must delete "Sites"
                   dataEnv =dune_trait_env$envir,
# delete "Species", "Species_abbr" from traits and
# use all remaining variables due to formulaTraits = ~. (the default)
                   dataTraits =dune_trait_env$traits[,-c(1,2)],
                   verbose = TRUE)
# Manual forward selection of environmental variables ----------------------------

# step 1
# user define:
consider <- names(dune_trait_env$envir)[c(2:6)]
test <- TRUE
cntr <- permute::how(nper= 999)
p.adjust.method <- "holm"
set.seed(213)


names(consider) <- consider
consider



# initiate
fit_measuresL <- list()
fit_measures <- matrix(NA, nrow = length(consider), ncol = 2)
colnames(fit_measures) <- c("eig1","pval1")
rownames(fit_measures) <- consider
considered <- NULL


for (k in seq_along(consider)){

  formulaE_FS <- as.formula(paste("~", consider[k]))

  out_FS <- dc_CA_vegan(formulaE_FS, dc_CA_vegan_object = out1, verbose = FALSE)

  if (test) {
    an <- anova(out_FS$RDAonEnv, permutations = cntr)
    pval <- an$`Pr(>F)`[1]} else pval <- NA
  fit_measures[k,] <- c(out_FS$eigenvalues[1], pval)
}

pvaladj <- p.adjust(fit_measures[,"pval1"], method = p.adjust.method)
fit_measures <- cbind(fit_measures,pvaladj )
round(fit_measures, 3) # best = Moist ; significant, also after correction for multiple testing....




chosen_env <- "Moist"


# step 2 ------------------------------------------------------------------
fit_measuresL[[chosen_env]] <- fit_measures
considered <- c(considered, consider[chosen_env])

considerk <- consider[-which(consider  %in% chosen_env)] # remove Moist


fit_measures <- matrix(NA, nrow = length(considerk), ncol = 2)
colnames(fit_measures) <- c("eig1","pval1")
rownames(fit_measures) <- considerk


for (k in seq_along(considerk)){

  formulaE_FS <- as.formula(paste("~",
      considerk[k], "+Condition(", paste(considered, collapse = "+") ,")", sep ="" ))

  out_FS <- dc_CA_vegan(formulaE_FS, dc_CA_vegan_object = out1, verbose = FALSE)

  if (test) {
    an <- anova(out_FS$RDAonEnv, permutations = cntr)
    pval <- an$`Pr(>F)`[1]} else pval <- NA
  fit_measures[k,] <- c(out_FS$eigenvalues[1], pval)
}

pvaladj <- p.adjust(fit_measures[,"pval1"], method = p.adjust.method)
fit_measures <- cbind(fit_measures,pvaladj )
round(fit_measures, 3) # best =Manure ; significant, also after correction for multiple testing (given step1)
chosen_env <- "Manure"

# step 3 ------------------------------------------------------------------
fit_measuresL[[chosen_env]] <- fit_measures
considered <- c(considered, consider[chosen_env])
considerk <- considerk[-which(considerk  %in% chosen_env)]

fit_measures <- matrix(NA, nrow = length(considerk), ncol = 2)
colnames(fit_measures) <- c("eig1","pval1")
rownames(fit_measures) <- considerk


for (k in seq_along(considerk)){
  formulaE_FS <- as.formula(paste("~", considerk[k],
      "+Condition(", paste(considered, collapse = "+") ,")", sep ="" ))

  out_FS <- dc_CA_vegan(formulaE_FS, dc_CA_vegan_object = out1, verbose = FALSE)

  if (test) {
    an <- anova(out_FS$RDAonEnv, permutations = cntr)
    pval <- an$`Pr(>F)`[1]} else pval <- NA
  fit_measures[k,] <- c(out_FS$eigenvalues[1], pval)
}

pvaladj <- p.adjust(fit_measures[,"pval1"], method = p.adjust.method)
fit_measures <- cbind(fit_measures,pvaladj )
round(fit_measures, 3) # best =Mag; not significant
chosen_env <- "Mag"

fit_measuresL[[chosen_env]] <- fit_measures


fit_measuresL
