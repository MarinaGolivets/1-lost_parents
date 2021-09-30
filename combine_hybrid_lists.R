# BiCIKL Hackathon (Sept 20-24, 2021, Meise)
# Combine all hybrid lists and match taxon names to the GBIF backbone


# load packages -------------------------------------------------------------------------------

pkgs <- c("here", "tidyverse", "parallel", "pbapply", "magrittr")
sapply(pkgs, require, character.only = TRUE)


# initialize parallelization ------------------------------------------------------------------

no_cores <- parallel::detectCores()
cl <- parallel::makeCluster(no_cores)


# load data sets ------------------------------------------------------------------------------

gbif_occ <- read_csv2(here("data/hybrids_gbif_occurrences.csv"))
gbif_chl <- read_csv2(here("data/hybrids_gbif_checklists.csv"))
wiki <- read_csv2(here("data/hybrids_wikidata.csv"))
ecoflora <- read_csv2(here("data/hybrids_ecoflora.csv"))
ipni <- read_csv2(here("data/hybrids_ipni.csv"))
floras <- read_csv2(here("data/hybrids_floras.csv"))
# vascan <- readLines(here("data/vascan_hybrid_names.txt"))
vascan <- read_csv2(here("data/vascan_hybrid_names.csv"), col_names = FALSE) %>%
  rename(verbatim = X1) %>%
  mutate(edited = gsub("×", " × ", verbatim)) %>%
  mutate(edited = str_trim(str_squish(edited)))

hybrids <- bind_rows(gbif_occ, gbif_chl, wiki, ecoflora, ipni, floras, vascan) %>%
  select(-id, -parent_taxonID) %>%
  distinct() %>%
  write_csv2(here("data/all_hybrids_raw.csv"))


# match to GBIF -------------------------------------------------------------------------------

source(here("fn_match_names_to_gbif.R"))
gbif.cmpfn <- compiler::cmpfun(gbif.fn)

matched <- gbif.cmpfn(
  taxonName = unique(hybrids$edited),
  taxonID = seq_along(unique(hybrids$edited))
) %>%
  select(taxonName, confidence, matchType, scientificName, usageKey, acceptedUsageKey) %>%
  mutate(acceptedUsageKey = ifelse(is.na(acceptedUsageKey), usageKey, acceptedUsageKey)) %>%
  rename(edited = taxonName) %>%
  distinct()

# # names with multiple IDs
# matched %>%
#   select(scientificName, usageKey) %>%
#   distinct() %>%
#   count(.$scientificName) %>%
#   filter(n > 1)

# exclude matches that are most likely incorrect
# View(matched %>%
#   filter((
#     grepl(" subsp\\.", scientificName) &
#       !grepl("subsp\\.|var\\.", edited) &
#       matchType == "FUZZY")))
# View(matched %>%
#   filter((
#     grepl(" var\\.", scientificName) &
#       !grepl("×", scientificName) &
#       !grepl("var\\.|subsp\\.", edited) &
#       matchType == "FUZZY")))

matched %<>%
  filter(!(
    grepl(" subsp\\.", scientificName) &
      !grepl("subsp\\.|var\\.", edited) &
      matchType == "FUZZY"))
matched %<>%
  filter(!(
    grepl(" var\\.", scientificName) &
      !grepl("×", scientificName) &
      !grepl("var\\.|subsp\\.", edited) &
      matchType == "FUZZY"))


# match parent taxa
parent_matched <- gbif.cmpfn(
  taxonName = unique(hybrids$parent),
  taxonID = seq_along(unique(hybrids$parent))
) %>%
  select(taxonName, confidence, matchType, scientificName, usageKey, acceptedUsageKey) %>%
  mutate(acceptedUsageKey = ifelse(is.na(acceptedUsageKey), usageKey, acceptedUsageKey)) %>%
  rename_with(~ paste0("parent_", .), .cols = everything()) %>%
  rename(parent = parent_taxonName) %>%
  distinct()

# # names with multiple IDs
# parent_matched %>%
#   select(parent_scientificName, parent_usageKey) %>%
#   distinct() %>%
#   count(.$parent_scientificName) %>%
#   filter(n > 1)

hybrids_matched <- hybrids %>%
  left_join(matched) %>%
  left_join(parent_matched) %>%
  filter(!is.na(usageKey)) %>%
  select(-verbatim, -parent) %>%
  distinct()

hybrids %<>%
  left_join(matched) %>%
  left_join(parent_matched) %>%
  filter(is.na(usageKey)) %>%
  bind_rows(hybrids_matched) %>%
  distinct()

w_parents <- hybrids %>%
  filter(!is.na(parent_scientificName))
wo_parents <- hybrids %>%
  filter(is.na(parent_scientificName)) %>%
  filter(!(scientificName %in% w_parents$scientificName))
hybrids <- bind_rows(w_parents, wo_parents)

n_distinct(hybrids$acceptedUsageKey) # 20999
n_distinct(hybrids %>% filter(is.na(acceptedUsageKey)) %>% select(edited)) # 25182

# matched hybrids with matched parents
hybrids %>% 
  filter(!is.na(acceptedUsageKey) & !is.na(parent_acceptedUsageKey)) %>%
  select(acceptedUsageKey, parent_acceptedUsageKey) %>%
  distinct() %>%
  pull(acceptedUsageKey) %>%
  n_distinct() # 6496

# non-matched hybrids with parents
hybrids %>% 
  filter(is.na(acceptedUsageKey) & !is.na(parent)) %>%
  select(edited, parent) %>%
  distinct() %>%
  pull(edited) %>%
  n_distinct() # 25182

write_csv2(hybrids, here("data/hybrids_master_list.csv"))