#' @title Dune meadow data with plant species traits and environmental variables

#' @description
#' The data \code{dune_trait_env} contains three data frames with abundance 
#' data of 28 plant species in 20 samples (relevés), trait data (9 traits:
#' of which 5 morphological and 4 ecological (Ellenberg indicator values) and
#' environmental data (9 environmental variables, four of which are geographic 
#' coordinates). Compared to the data in Jongman et al. (1987, 1995), the two 
#' moss species are lacking, and the traits of plant species and the geographic 
#' coordinates of the samples are added. The data and the following description 
#' are an edited version of the DataKey in the Jamil2013_AJ data set in the 
#' CESTES database (Jeliazkov et al. 2020).
#'
#' The Dune Meadow Data originate from a MSc thesis report of Batterink & 
#' Wijffels (1983). It consisted of 80 relevés in about 68 dune meadows (lots) 
#' on the island Terschelling in the Netherlands. A subset of their data was 
#' selected by Caspar Looman as an example data set for the edited book 
#' Jongman, ter Braak and van Tongeren (1987, 1995). The subset consists of 20 
#' relevés, 28 species (and 2 mosses, excluded here) and 5 environmental 
#' variables.
#' 
#' The trait data were taken from the LEDA database by  Jamil et al 2013. The 
#' spatial coordinates were retrieved by geo-referencing in GIS of the maps in
#' Batterink & Wijffels (1993) by Ruut Wegman and Cajo ter Braak. The X, Y 
#' coordinates are by geo-referencing the relevé locations on Kaart 3a, 3b and 
#' 3c; the X_lot, Y_lot coordinates are from Kaart 2a, 2b and 2c.
#' 
#' For the Ellenberg indicator values see Ellenberg (1992).
#'
#' The data \code{dune_trait_env} is a list with elements that are data frames 
#' each
#' \itemize{
#' \item \code{comm}:   community data; vegetation data.
#' \item \code{traits}: trait data, taken from the LEDA database.
#' \item \code{envir}:  environmental data, taken from Jongman et al. (1987,1995).
#' }
#' The community data collection was done by the Braun-Blanquet method;
#' the data are recorded according to the ordinal scale of van der Maarel (1979, 
#' Vegetatio, 39, 97-114); see pages XVII-XVIII and 18 in Jongman, ter Braak & 
#' van Tongeren 1995. Nomenclature follows Heukels-Van der Meijden (1983) 
#' Flora van Nederland, 20th ed.
#'
#' The morphological \code{traits} are
#' \itemize{
#' \item \code{SLA}:	    Specific Leaf Area
#' \item \code{Height}:	  Canopy height of a shoot
#' \item \code{LDMC}:  	  Leaf dry matter content
#' \item \code{Seedmass}:	Seed mass
#' \item \code{Lifespan}:	Life span. Nominal; annual vs. perennial
#' }
#' The ecological \code{traits} (habitat requirements) are the Ellenberg values
#' \itemize{
#' \item \code{F}:	    Moisture (ranging [1 to 12] (low to high))
#' \item \code{R}:	    Soil acidity, ranging [1 to 9] (acidic to alkaline)
#' \item \code{N}:  	  Nitrogen requirement, ranging [1 to 9] (low to high)
#' \item \code{L}:	    Light requirement, ranging [1 to 9] (low to high)
#' } 
#'  
#' The data frame \code{envir} contains the environmental variables
#' \itemize{
#' \item \code{A1}:    horizon thickness
#' \item \code{Moist}:	Moisture content of the soil (a five point scale)
#' \item \code{Mag}:    Grassland management type
#' \item \code{Use}:	type of use (Agricultural grassland use (1) hay production 
#' (2) intermediate (3) grazing)
#' \item \code{Manure}:	Quantity of manure applied based on N and P manuring 
#' (N/P class in B&W 1983)
#' \item \code{X}:	longitude	geographical coordinates (m) of the 2x2 m2 
#' sample (relevé)
#' \item \code{Y}:	latitude	geographical coordinates (m) of the 2x2 m2 
#' sample (relevé)
#' \item \code{X_lot}:	longitude	geographical coordinates (m) of the lot center
#' \item \code{Y_lot}:  latitude	geographical coordinates (m) of the lot center
#' }
#' The management types are standard farming (SF), biological farming (BF),
#' hobby farming (HF), nature conservation management (NM). The coordinates are
#' Rijksdriehoekscoordinaten in meters.
#' \url{https://nl.wikipedia.org/wiki/Rijksdriehoekscoordinaten}
#' 
#' @name dune_trait_env
#' @docType data
#' @references
#' Batterink, M. & Wijffels, G. (1983) Een vergelijkend vegetatiekundig 
#' onderzoek naar de typologie en invloeden van het beheer van 1973 tot 1982 in 
#' de duinweilanden op Terschelling. Landbouwhogeschool. 
#' ISN 215909.01. WUR library stacks 704B58.
#' 
#' Ellenberg, H. (1992) Indicator values of plants in Central Europe.
#' Scripta Geobotanica, 18, 258 pp.
#'
#' Jamil, T., Ozinga, W.A., Kleyer, M. & ter Braak, C.J.F. (2013) Selecting 
#' traits that explain species–environment relationships: a generalized linear 
#' mixed model approach. Journal of Vegetation Science, 24, 988-1000.
#' \doi{10.1111/j.1654-1103.2012.12036.x}.
#'
#' Jeliazkov, A., Mijatovic, D., and 78 others. (2020) A global database for 
#' metacommunity ecology, integrating species, traits, environment and space. 
#' Scientific Data, 7. \doi{10.1038/s41597-019-0344-7}.
#'
#' Jongman, R.H.G., ter Braak, C.J.F. & van Tongeren, O.F.R. (1987) Data 
#' analysis in community and landscape ecology. Pudoc, Wageningen. 
#' ISBN 90-220-0908-4.
#'
#' Jongman, R.H.G., ter Braak, C.J.F. & van Tongeren, O.F.R. (1995) Data 
#' analysis in community and landscape ecology. Cambridge University Press, 
#' Cambridge. ISBN 0-521-47574-0. \url{https://edepot.wur.nl/248017}.
#'
#' Kleyer, M., and 33 others (2008) The LEDA Traitbase: a database of 
#' life-history traits of the Northwest European flora. Journal of Ecology, 
#' 96, 1266-1274. \doi{10.1111/j.1365-2745.2008.01430.x}.
NULL
