#' @title Extract results of a wrda or dcca object
#'
#' @description
#' This function works very much like the \code{vegan} 
#' \code{\link[vegan]{scores}} function, in particular 
#' \code{\link[vegan]{scores.cca}}, with the additional results such as 
#' regression coefficients and linear combinations of traits 
#' \code{('regr_traits','lc_traits')}. All scores from CA obey the so called
#' transition formulas and so do the scores of CCA and dc-CA. The differences
#' are, for CCA, that the linear combinations of environmental variables (the 
#' \emph{constrained} site scores) replace the usual (\emph{unconstrained}) 
#' site scores, and for dc-CA, that the linear combinations of traits (the 
#' \emph{constrained} species scores) also replace the usual 
#' (\emph{unconstrained}) species scores.
#'
#' @param x object of class \code{"dcca"}, \emph{i.e.} result of
#' \code{\link{dc_CA}}.
#' @param choices integer vector of which axes to obtain. Default: all dc-CA 
#' axes.
#' @param display a character vector, one or more of \code{c("all", "species",
#' "sites", "sp", "wa", "lc", "bp", "cor", "ic", "reg", "tval", "cn", 
#' "lc_traits", "reg_traits","tval_traits", "cor_traits", "ic_traits", 
#' "bp_traits","cn_traits")}. The most items are as in 
#' \code{\link[vegan]{scores.cca}}, except \code{"cor"} and \code{"ic"}, for
#' inter-set and intra-set correlations, respectively, and \code{"tval"} for
#' the (over-optimistic) t-values of the regression coefficients. The remaining
#' scores are analogous scores for species and traits.
#' @param which_cor character or list of trait and environmental variables 
#' names (in this order) in the data frames for which inter-set correlations 
#' must calculated. Default: a character ("in_model") for all traits and 
#' variables in the model, including collinear variables and levels.
#' @param scaling numeric (1,2 or 3) or character \code{"sites", "species" or
#' "symmetric"}. Default: "symmetric". Either site- (1) or species- (2) related
#' scores are scaled by eigenvalues, and the other set of scores is left 
#' unscaled, or with 3 both are scaled symmetrically by square root of 
#' eigenvalues. Negative values are treated as the corresponding positive ones
#' by \code{abs(scaling)}.
#' @param tidy Return scores that are compatible with \code{ggplot2}: all 
#' scores are in a single data.frame, score type is identified by factor 
#' variable \code{score}, the names by variable \code{label}, and species 
#' weights (in \code{\link{dc_CA}} are in variable \code{weight}. See 
#' \code{\link[vegan]{scores.cca}}.
#' @param ... Other arguments passed to the function (currently ignored).
#' 
#' @details
#' The function is modeled after \code{\link[vegan]{scores.cca}}.
#'
#' If you get the error message: 'arg' should be one of "sites", "species", 
#' "both", then the vegan scores function has been called, instead of the one 
#' of douconca. The work-around is to use douconca::scores() instead of 
#' scores() only.
#'
#' An example of which_cor is: \code{which_cor = list(traits = "SLA", 
#' env = c("acidity", "humidity"))}
#' 
#' @return A data frame if \code{tidy = TRUE}, a matrix if a single item is 
#' asked for and a named list of matrices if more than one item is asked for. 
#' The following names can be included: \code{c("sites", "constraints_sites",
#' "centroids", "regression", "t_values", "correlation", 
#' "intra_set_correlation", "biplot", "species", "constraints_species", 
#' "regression_traits", "t_values_traits", "correlation_traits",
#' "intra_set_correlation_traits", "biplot_traits", "centroids_traits")}. Each
#' matrix has an attribute \code{"meaning"} explaining its meaning. With 
#' \code{tidy = TRUE}, the resulting data frame has attributes \code{"scaling"} 
#' and \code{"meaning"}; the latter has two columns: (1) name of score type 
#' and (2) its meaning, usage and interpretation.
#'
#' An example of the meaning of scores in scaling \code{"symmetric"} with 
#' \code{display ="all"}:
#' \describe{
#'  \item{sites}{ CMWs of the trait axes (constraints species) in scaling 
#'  'symmetric' optimal for biplots and, almost so, for inter-site distances.}
#'  \item{constraints_sites}{linear combination of the environmental predictors 
#'  and the covariates (making the ordination axes orthogonal to the 
#'  covariates) in scaling 'symmetric' optimal for biplots and, almost so, 
#'  for inter-site distances.}
#'  \item{regression}{mean, sd, VIF, standardized regression coefficients and 
#'  their optimistic t-ratio in scaling 'symmetric'.}
#'  \item{t_values}{t-values of the coefficients of the regression of the CWMs 
#'  of the trait composite on to the environmental variables}
#'  \item{correlation}{inter set correlation, correlation between environmental 
#'  variables and the sites scores (CWMs)}
#'  \item{intra_set_correlation}{intra set correlation, correlation between
#'  environmental variables and the dc-ca axis (constrained sites scores)}
#'  \item{biplot}{biplot scores of environmental variables for display with 
#'  biplot-traits for fourth-corner correlations in scaling 'symmetric'.}
#'  \item{centroids}{environmental category means of the site scores in scaling 
#'  'symmetric'  optimal for biplots and, almost so, for inter-environmental
#'  category distances.}
#'  \item{species}{SNC on the environmental axes (constraints sites) in scaling 
#'  'symmetric' optimal for biplots and, almost so, for inter-species 
#'  distances.}
#'  \item{constraints_species}{linear combination of the traits and the trait 
#'  covariates (making the ordination axes orthogonal to the covariates) in 
#'  scaling 'symmetric' optimal for biplots and, almost so, for inter-species 
#'  distances.}
#'  \item{regression_traits}{mean, sd, VIF, standardized regression 
#'  coefficients and their optimistic t-ratio in scaling 'symmetric'.}
#'  \item{t_values_traits}{t-values of the coefficients of the regression of the
#'  SNCs along a dc-CA axis on to the traits}
#'  \item{correlation_traits}{inter set correlation, correlation between 
#'  traits and the species scores (SNCs)}
#'  \item{intra_set_correlation_traits}{intra set correlation, correlation 
#'  between traits and the dc-ca axis (constrained species scores)}
#'  \item{biplot_traits}{biplot scores of traits for display with biplot scores
#'  for fourth-corner correlation in scaling 'symmetric'.}
#'  \item{centroids_traits}{trait category means of the species scores in 
#'  scaling 'symmetric' optimal for biplots and, almost so, for inter-trait 
#'  category distances.}
#' }
#'
#' The statements on optimality for distance interpretations are based on the 
#' \code{scaling} and the relative magnitude of the dc-CA eigenvalues of the 
#' chosen axes.
#' 
#' @example demo/dune_dcCA.R
#' 
#' @noRd
#' @keywords internal
scores_dcca <- function(x, 
                        choices = 1:2, 
                        display = "all", 
                        scaling = "sym", 
                        which_cor = "in model", 
                        tidy = FALSE,
                        ...) {
  # internal function
  f_meaning <- function(type_of_scores, scaling, txt) {
    if (type_of_scores %in% c("sites", "constraints")) {
      point_type <- "site" 
      opt_scal <- "sites"
    } else if (type_of_scores %in% c("species", "constraints_species")) {
      point_type <- "species" 
      opt_scal <- "species"
    } else if (type_of_scores %in% "centroids") {
      point_type <- "environmental category" 
      opt_scal <- "sites"
    } else if (type_of_scores %in% "centroids_traits") {
      point_type <- "trait category"
      opt_scal <- "species"
    } else if (type_of_scores %in% c("biplot", "biplot_traits")) {
      point_type <- "arrows"
      opt_scal <- scaling
    }
    txt1 <- paste0(txt, " in scaling '", scaling, 
                   "' optimal for biplots and inter-", point_type, " distances.")
    txt1b <- paste0(txt, " in scaling '", scaling, 
                    "' optimal for biplots and, almost so, for inter-", 
                    point_type, " distances.")
    txt1c <- paste0(txt, " in scaling '", scaling, 
                    "' optimal for biplots, but unsuited for inter-", 
                    point_type, " distances.")
    txt2 <- paste0(txt, " in scaling '", scaling, 
                   "', optimal for biplots displays, suboptimal for display of inter-",
                   point_type, " distances.")
    txt3 <- paste0(txt, " in scaling '", scaling, 
                   "' optimal for biplot displays, suboptimal for distance interpretation.")
    txt4 <- paste0(txt, " in scaling '", scaling, "'.")
    if (type_of_scores %in% c("biplot", "biplot_traits")) {
      txt1 <- txt2 <- txt3 <- txt4
    }
    thres1 <- 1.5 # to be decided upon... What does Canoco 5 suggest?
    thres2 <- 4
    if (length(choices) > 1) {
      ratio_eig <- x$eigenvalues[choices[1]] / x$eigenvalues[choices[2]]
      if (scaling == "symmetric") {
        if (sqrt(ratio_eig) < thres1) {
          txt3 <- txt1b 
        } else if (sqrt(ratio_eig) > thres2) {
          txt3 <- txt1c
        }
      } else {
        if (ratio_eig < thres1) {
          txt2 <- txt1b 
        } else if (ratio_eig > thres2) {
          txt2 <- txt1c
        }
      }
      if (opt_scal == scaling) {
        txt_out <- txt1
      } else if (scaling == "symmetric") {
        txt_out <- txt3
      } else {
        txt_out <- txt2
      }
    } else txt_out <- txt
    return(txt_out)
  }
  tabula <- c("sites", "constraints", "regression", "t_values", "biplot", 
              "correlation", "intra_set_correlation", "centroids", "species",
              "constraints_species", "regression_traits", "t_values_traits",
              "biplot_traits", "correlation_traits", 
              "intra_set_correlation_traits", "centroids_traits")
  names(tabula) <- c("wa", "lc", "reg", "tval", "bp", "cor", "ic", "cn", "sp",
                     "lc_traits", "reg_traits", "tval_traits", "bp_traits",
                     "cor_traits", "ic_traits", "cn_traits")
  arg <- c("sp", "wa", "lc", "bp", "cor", "ic", "reg", "tval", "cn", 
           "lc_traits", "reg_traits", "tval_traits", "bp_traits", 
           "cor_traits", "ic_traits", "cn_traits", "sites", "species", "all")
  if (tidy) {
    display <- "all"
  }
  if (inherits(x, "dcca")) {
    display <- match.arg(display, arg, several.ok = TRUE)
    if ("all" %in% display) {
      display <- names(tabula)
    }
  } else if (inherits(x, "wrda", which = TRUE) == 1) {
    display <- match.arg(display, arg[-c(10:16)], several.ok = TRUE)
    if ("all" %in% display) {
      display <- names(tabula)[-c(10:16)]
    } 
    if (!which_cor[1] == "in model") {
      which_cor <- list(NA, which_cor)
    }
  } else {
    stop("The first argument must be of class 'dcca' or 'wrda', the result ", 
         "of the function dc_CA or wrda.")
  }
  ## set "all" for tidy scores
  if ("sites" %in% display) {
    display[display == "sites"] <- "wa"
  }
  if ("species" %in% display) {
    display[display == "species"] <- "sp"
  }
  if ("correlation" %in% display) {
    display[display == "correlation"] <- "co"
  }
  take <- tabula[display]
  if (inherits(x, "wrda", which = TRUE) == 1) {
    site_axes <- f_env_axes(x)
    species_axes <- x$species_axes
    formulaEnv <- x$formula
    dataEnv <- x$data
  } else { # dcca
    dataEnv <- x$data$dataEnv
    formulaEnv <- x$formulaEnv
    if (!"species_axes" %in% names(x)) {
      site_axes <- f_env_axes(x)
      species_axes <- f_trait_axes(x)
    } else {
      species_axes <- x$species_axes
      site_axes<- x$site_axes
    }
  }
  # make sure axes chosen by choices are not larger than the rank
  choices <- choices[choices <= length(x$eigenvalues)]
  if (is.character (scaling)) {
    scaling <- match.arg(scaling,  c("sites", "species", "symmetric"))
    if (scaling == "sites") {
      num_scaling <- 1 
    } else if (scaling == "species") {
      num_scaling <- 2 
    } else if (scaling == "symmetric") {
      num_scaling <- 3 
    }
  } else if (is.numeric(scaling)) {
    num_scaling <- abs(scaling)
    scaling <- c("sites","species","symmetric")[num_scaling]
  } else {
    stop("scaling type not recognized")
  }
  slam <- sqrt(x$eigenvalues[choices])
  scal <- list(rep(1, length(slam)), slam, sqrt(slam))[[abs(num_scaling)]]
  if (length(scal) > 1) {
    diag_scal_sites <- diag(1 / scal)
    diag_scal_species <- diag(scal)
  } else {
    diag_scal_sites   <- matrix(1 / scal)
    diag_scal_species <- matrix(scal)
  }
  sol <- list()
  # scaling for site related scores (incl env)
  if ("sites" %in% take) {
    sol$sites <- 
      site_axes$site_scores[[1]][, choices, drop = FALSE] %*% diag_scal_sites
    attr(sol$sites, which = "meaning") <- 
      f_meaning("sites", scaling, 
                "CMWs of the trait axes (constraints species)")
  }
  if ("constraints" %in% take) {
    sol$constraints_sites <- 
      site_axes$site_scores[[2]][,choices, drop = FALSE] %*% diag_scal_sites
    attr(sol$constraints_sites, which = "meaning") <- 
      f_meaning("constraints", scaling,
                paste("linear combination of the environmental predictors",
                      "and the covariates (making the ordination axes", 
                      "orthogonal to the covariates)"))
  }
  if ("regression" %in% take) {
    regr <- 
      site_axes$c_env_normed[, choices + 3, drop = FALSE] %*% diag_scal_sites
    if (tidy) {
      sol$regression <- regr 
    } else {
      sol$regression <- cbind(site_axes$c_env_normed[, 1:3, drop = FALSE], regr)
    }
    attr(sol$regression, which = "meaning")<-
      paste0("mean, sd, VIF, standardized regression coefficients and their ", 
             "optimistic t-ratio in scaling '", scaling, "'.")
  }
  if ("t_values" %in% take){
    sol$t_values <- 
      site_axes$c_env_normed[, Rank_mod(x) + choices + 3, drop = FALSE]
    attr(sol$t_values, which = "meaning") <-
      paste("t-values of the coefficients of the regression of the CWMs", 
            "of the trait composite on to the environmental variables")
  }
  if ("correlation" %in% take) {
    if (!is.list(which_cor)) {
      sol$correlation <- site_axes$correlation[, choices, drop = FALSE]
    } else {
      whichc <- which_cor[[2]]
      if (whichc[1] == "in model") {
        whichc <- get_Z_X_XZ_formula(formulaEnv)$focal_nams
      }
      cor_Env_CWM <- f_env_axes(x, which_cor = whichc)
      sol$correlation <- cor_Env_CWM$correlation[, choices, drop = FALSE]
    }
    attr(sol$correlation, which = "meaning") <-
      paste("inter set correlation, correlation between environmental", 
            "variables and the sites scores (CWMs)")
  }
  if ("intra_set_correlation" %in% take) {
    if (!is.list(which_cor)) {
      e_rcor <- site_axes$correlation[, choices, drop = FALSE]
    } else {
      whichc <- which_cor[[2]]
      if (whichc[1] == "in model") {
        whichc <- get_Z_X_XZ_formula(formulaEnv)$focal_nams
      }
      cor_Env_CWM <- f_env_axes(x, which_cor = whichc)
      e_rcor <- cor_Env_CWM$correlation[, choices, drop = FALSE]
    }
    R <- sqrt(site_axes$R2_env[choices])
    if (length(R) == 1) {
      sR <- matrix(1 / R)
    } else {
      sR <- diag(1 / R)
    }
    sol$intra_set_correlation <- e_rcor %*% sR
    colnames(sol$intra_set_correlation) <- paste0("dcCA", choices)
    attr(sol$intra_set_correlation,  which = "meaning") <-
      paste("intra set correlation, correlation between environmental", 
            "variables and the dc-ca axis (constrained sites scores)")
  }
  if ("biplot" %in% take) {
    e_rcor <- site_axes$correlation[, choices, drop = FALSE]
    R <- sqrt(site_axes$R2_env[choices])
    if (length(R) == 1) {
      sR <- matrix(slam / R)
    } else {
      sR <- diag(slam / R)
    }
    sol$biplot <- e_rcor %*% sR %*% diag_scal_sites
    colnames(sol$biplot) <- paste0("dcCA", choices)
    attr(sol$biplot, which = "meaning") <- 
      f_meaning("biplot", scaling,
                paste("biplot scores of environmental variables for display", 
                      "with biplot-traits for fourth-corner correlations"))
  }
  if ("centroids" %in% take) {
    if (!is.list(which_cor)) {
      whichc <- get_Z_X_XZ_formula(formulaEnv, dataEnv)$focal_nams
    } else {
      whichc <- which_cor[[2]]
      if (whichc[1] ==  "in model") {
        whichc <- get_Z_X_XZ_formula(formulaEnv)$focal_nams
      }
    }
    dat <- dataEnv[, whichc, drop = FALSE]
    cn <- centroids.cca(site_axes$site_scores$site_scores_unconstrained,
                        dat, wt = x$weights$rows)[, choices, drop = FALSE]
    if (!is.null(cn)) {
      cn <- cn %*% diag_scal_sites
      attr(cn, which = "meaning") <- 
        f_meaning("centroids", scaling, 
                  "environmental category means of the site scores")
    }
    sol$centroids <-  cn
  }
  # Species stats
  if ("species" %in% take) {
    sol$species <- 
      species_axes$species_scores[[1]][, choices, drop = FALSE] %*% diag_scal_species
    attr(sol$species, which = "meaning")<-
      f_meaning("species", scaling, 
                "SNC on the environmental axes (constraints sites)")
  }
  if (inherits(x, "wrda", which = TRUE) != 1){
    # dcca
    if ("constraints_species" %in% take) {
      sol$constraints_species <- 
        species_axes$species_scores[[2]][, choices, drop = FALSE] %*% diag_scal_species
      attr(sol$constraints_species, which = "meaning") <- 
        f_meaning("constraints_species", scaling,
                  paste("linear combination of the traits and the trait",
                        "covariates (making the ordination axes orthogonal", 
                        "to the covariates)'")
        )
    }
    if ("regression_traits" %in% take) {
      regr <- 
        species_axes$c_traits_normed[, choices + 3, drop = FALSE] %*% diag_scal_species
      if (tidy) {
        sol$regression_traits <- regr 
      } else {
        sol$regression_traits <- 
          cbind(species_axes$c_traits_normed[, 1:3, drop = FALSE], regr)  
      }
      attr(sol$regression_traits, which = "meaning") <-
        paste0("mean, sd, VIF, standardized regression coefficients and their ", 
               "optimistic t-ratio in scaling '", scaling, "'.")
    }
    if ("t_values_traits" %in% take) {
      sol$t_values_traits <- 
        species_axes$c_traits_normed[, Rank_mod(x) + choices + 3, drop = FALSE]
      attr(sol$t_values_traits, which = "meaning") <-
        paste("t-values of the coefficients of the regression of the SNCs", 
              "along a dc-CA axis on to the traits")
    }
    if ("correlation_traits" %in% take) {
      if (!is.list(which_cor)) {
        sol$correlation_traits <- 
          species_axes$correlation[, choices, drop = FALSE]
      } else {
        whichc <- which_cor[[1]]
        if (whichc[1] == "in model") {
          whichc <- get_Z_X_XZ_formula(x$formulaTraits)$focal_nams
        }
        Cor_Trait_SNC <- f_trait_axes(x, which_cor = whichc)
        sol$correlation_traits <- 
          Cor_Trait_SNC$correlation[, choices, drop = FALSE]
      }
      attr(sol$correlation_traits, which = "meaning") <-
        paste("inter set correlation, correlation between traits and the",
              "species scores (SNCs)")
    }
    if ("intra_set_correlation_traits" %in% take) {
      if (!is.list(which_cor)) {
        e_rcor <- species_axes$correlation[, choices, drop = FALSE]
      } else {
        whichc <- which_cor[[1]]
        if (whichc[1] == "in model") {
          whichc <- get_Z_X_XZ_formula(x$formulaTraits)$focal_nams
        }
        Cor_Trait_SNC <- f_trait_axes(x, which_cor = whichc)
        e_rcor <- Cor_Trait_SNC$correlation[, choices, drop = FALSE]
      }
      R <- sqrt(species_axes$R2_traits[choices])
      if (length(R) == 1) {
        sR <- matrix(1 / R)
      } else {
        sR <- diag(1 / R)
      }
      sol$intra_set_correlation_traits <- e_rcor %*% sR
      colnames(sol$intra_set_correlation_traits) <- 
        paste0("dcCA", choices)
      attr(sol$intra_set_correlation_traits, which = "meaning") <-
        paste("intra set correlation, correlation between traits and the", 
              "dc-ca axis (constrained species scores)")
    }
    if ("biplot_traits" %in% take) {
      t_rcor <- species_axes$correlation[, choices, drop = FALSE]
      R <- sqrt(species_axes$R2_traits[choices])
      if (length(R) == 1) {
        sR <- diag_scal_species / R 
      } else {
        sR <-  diag(diag(diag_scal_species) / R)
      }
      sol$biplot_traits <- t_rcor %*% sR
      colnames(sol$biplot_traits) <- paste0("dcCA", choices)
      attr(sol$biplot_traits, which = "meaning") <-
        f_meaning("biplot", scaling, 
                  paste("biplot scores of traits for display with biplot", 
                        "scores for fourth-corner correlation"))
    }
    if ("centroids_traits" %in%take) {
      if (!is.list(which_cor)) {
        whichc <- get_Z_X_XZ_formula(x$formulaTraits)$focal_nams
      } else {
        whichc <- which_cor[[1]]
        if (whichc[1] == "in model") {
          whichc <- get_Z_X_XZ_formula(x$formulaTraits)$focal_nams
        }
      }
      dat <- x$data$dataTraits[, whichc, drop = FALSE]
      cn <- centroids.cca(species_axes$species_scores$species_scores_unconstrained,
                          dat, wt = x$weights$columns)[, choices, drop = FALSE]
      if (!is.null(cn)) {
        cn <- cn %*% diag_scal_species
        attr(cn, which = "meaning") <-
          f_meaning("centroids_traits", scaling,
                    "trait category means of the species scores")
      }
      sol$centroids_traits <-  cn
    }
  }
  # end of types of scores
  for (nam in names(sol)){
    if(!is.null(sol[[nam]])) {
      if (!nam %in% c("regression", "regression_traits", "correlation", 
                      "correlation_traits")) {
        colnames(sol[[nam]]) <- paste0("dcCA", choices) 
      } else if (nam %in% c("regression", "regression_traits")) {
        colnames(sol[[nam]])[-c(1, 2, 3)] <- paste0("dcCA", choices)
      }
    }
  }
  # taken from vegan::scores.cca with thanks
  ## Take care that scores have names
  if (length(sol)) {
    for (i in seq_along(sol)) {
      if (is.matrix(sol[[i]])) {
        rownames(sol[[i]]) <- rownames(sol[[i]], do.NULL = FALSE,
                                       prefix = substr(names(sol)[i], 1, 3))
      }
    }
  }
  # meaning of scores for tidy
  meaning0 <- lapply(sol, FUN = attr, which = "meaning")
  meaning1 <- character(length(meaning0))
  for (i in seq_along(meaning0)) {
    meaning1[i] <- meaning0[[i]]
  }
  meaning <- data.frame(type = names(meaning0), meaning = meaning1)
  ## tidy scores
  if (tidy) {
    if (length(sol) == 0) {# no requested scores existed
      return(NULL)
    }
    ## re-group biplot arrays duplicating factor centroids
    if (!is.null(sol$biplot) && !is.null(sol$centroids)) {
      dup <- rownames(sol$biplot) %in% rownames(sol$centroids)
      if (any(dup)) {
        sol$factorbiplot <- sol$biplot[dup, , drop = FALSE]
        sol$biplot <- sol$biplot[!dup, , drop = FALSE]
      }
    }
    ## re-group biplot arrays duplicating factor centroids
    if (!is.null(sol$biplot_traits) && !is.null(sol$centroids_traits)) {
      dup <- rownames(sol$biplot_traits) %in% rownames(sol$centroids_traits)
      if (any(dup)) {
        sol$factorbiplot_traits <- sol$biplot_traits[dup, , drop = FALSE]
        sol$biplot_traits <- sol$biplot_traits[!dup, , drop = FALSE]
      }
    }
    group <- sapply(sol, nrow)
    group <- rep(names(group), group)
    sol <- do.call(rbind, sol)
    label <- rownames(sol)
    cw <- x$weights$columns
    rw <- x$weights$rows
    w <- rep(NA, nrow(sol))
    if (any(weighted <- group == "sites")) {
      w[weighted] <- rw
    }
    if (any(weighted <- group == "constraints")) {
      w[weighted] <- rw
    }
    if (any(weighted <- group == "species")) {
      w[weighted] <- cw
    }
    if (any(weighted <- group == "constraints_species")) {
      w[weighted] <- cw
    }
    sol <- as.data.frame(sol)
    sol$score <- as.factor(group)
    sol$label <- label
    sol$weight <- w
    names(sol)[seq_along(choices)] <- paste0("dcCA", choices)
    attr(sol, which = "scaling") <- scaling
    attr(sol, which = "meaning") <- meaning
  }
  ## return NULL instead of list(), and matrix instead of a list of
  ## one matrix
  return(switch(min(2, length(sol)), sol[[1]], sol))
}
