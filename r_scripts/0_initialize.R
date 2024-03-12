# Marina Golivets
# BiCIKL Hackathon (Sept 20-24, 2021, Meise, Belgium)
# Initialize 


# load packages ------------------------------------------------------------------------------------

install.packages("remotes")
remotes::install_github("barnabywalker/kewr")

pkgs <- c("here", "tidyverse", "magrittr", "janitor", "parallel", 
          "data.table", "foreach", "stringi", "WorldFlora")
sapply(pkgs, require, character.only = TRUE)


# initialize parallelization -----------------------------------------------------------------------
no_cores <- parallel::detectCores()
cl <- parallel::makeCluster(no_cores)

source(here("fn_match_names_to_gbif.R"))
gbif.cmpfn <- compiler::cmpfun(gbif.fn)
source(here("fn_gbif_get-accepted-names.R"))


parallel::stopCluster(cl)

