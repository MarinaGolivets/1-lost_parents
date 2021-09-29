# Marina Golivets
# BiCIKL Hackathon (Sept 20-24, 2021, Meise)
# Extract hybrid and parent taxon names from Wikidata


# load packages ------------------------------------------------------------------------------

pkgs <- c("here", "tidyverse", "magrittr")
sapply(pkgs, require, character.only = TRUE)


# initialize parallelization ------------------------------------------------------------------

no_cores <- parallel::detectCores()
cl <- parallel::makeCluster(no_cores)


# load taxon names from Wikidata --------------------------------------------------------------

hybrids_wiki <- readxl::read_excel(here("data/Summary_Wikidata.xlsx"), sheet = 2, skip = 1) %>%
  janitor::clean_names() %>%
  rename(verbatim = item_label_2, parent = item_label_4) %>%
  select(verbatim, parent)


# match to GBIF backbone ---------------------------------------------------------------------

source(here("fn_match_names_to_gbif.R"))
gbif.cmpfn <- compiler::cmpfun(gbif.fn)

matched <- gbif.cmpfn(
  taxonName = unique(hybrids_wiki$verbatim),
  taxonID = seq_along(unique(hybrids_wiki$verbatim))
) %>%
  select(taxonName, scientificName, matchType, confidence, usageKey) %>%
  rename(verbatim = taxonName)

parent_matched <- gbif.cmpfn(
  taxonName = unique(hybrids_wiki$parent),
  taxonID = seq_along(unique(hybrids_wiki$parent))
) %>%
  select(taxonName, scientificName, usageKey, matchType, confidence) %>%
  rename_with(~ paste0("parent_", .), .cols = everything()) %>%
  mutate(parent = parent_taxonName)

hybrids_wiki %<>%
  left_join(matched) %>%
  left_join(parent_matched) %>%
  filter(!is.na(scientificName)) %>%
  mutate(parent = ifelse(is.na(parent_taxonName), NA, parent)) %>%
  select(verbatim, parent)

w_parents <- hybrids_wiki %>%
  filter(!is.na(parent))
wo_parents <- hybrids_wiki %>%
  filter(is.na(parent)) %>%
  filter(!(verbatim %in% w_parents$verbatim))
hybrids_wiki <- bind_rows(w_parents, wo_parents)

write_csv2(hybrids_wiki, here("data/hybrids_wikidata.csv"))