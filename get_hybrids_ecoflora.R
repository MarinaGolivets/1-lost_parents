# Marina Golivets
# BiCIKL Hackathon (Sept 20-24, 2021, Meise)
# Extract hybrid and parent taxon names from EcoFlora


# load packages ------------------------------------------------------------------------------------

pkgs <- c("here", "tidyverse", "magrittr")
sapply(pkgs, require, character.only = TRUE)


# get hybrid names from EcoFlora ------------------------------------------------------------------

hybrids_ecoflora <- read_csv2(here("data/ecoflora_taxa.csv")) %>%
  filter(is_hybrid == TRUE) %>%
  filter(!grepl(" '", taxon_name)) %>%
  rename(verbatim = taxon_name, formula = hybrid_formula) %>%
  mutate(formula = ifelse(is.na(formula), verbatim, formula)) %>%
  mutate(parent = str_split(formula, pattern = " x ")) %>%
  unnest(cols = "parent") %>%
  rowwise() %>%
  mutate(genus = str_split(verbatim, pattern = " ", simplify = TRUE)[1]) %>%
  mutate(parent = gsub("^[A-Z]\\.", genus, parent)) %>%
  mutate(parent = ifelse(
    grepl("^[a-z]", parent),
    paste(genus, parent), parent
  )) %>%
  select(verbatim, parent) %>%
  mutate(edited = gsub(" x |^X ", " × ", verbatim)) %>%
  mutate(edited = str_trim(str_squish(edited))) %>%
  ungroup() %>%
  mutate(parent = ifelse(!grepl(" ", parent), NA, parent)) %>%
  distinct() %>%
  write_csv2(here("data/hybrids_ecoflora.csv"))
