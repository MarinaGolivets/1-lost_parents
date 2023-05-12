# Marina Golivets
# BiCIKL Hackathon (Sept 20-24, 2021, Meise, Belgium)
# Initialize 


# load packages ------------------------------------------------------------------------------------
pkgs <- c("here", "tidyverse", "magrittr")
sapply(pkgs, require, character.only = TRUE)

# initialize parallelization -----------------------------------------------------------------------
no_cores <- parallel::detectCores()
cl <- parallel::makeCluster(no_cores)

source(here("fn_match_names_to_gbif.R"))
gbif.cmpfn <- compiler::cmpfun(gbif.fn)

parallel::stopCluster(cl)
parallel::stopCluster(cl)

