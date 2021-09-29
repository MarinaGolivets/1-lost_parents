# Marina Golivets
# BiCIKL Hackathon (Sept 20-24, 2021, Meise)
# Extract hybrid and parent taxon names from EcoFlora


# load packages ------------------------------------------------------------------------------------

pkgs <- c("here", "tidyverse", "magrittr")
sapply(pkgs, require, character.only = TRUE)


# get hybrid names from IPNI -------------------------------------------------------------------

hybrids_ipni <- data.table::fread(here("data/ipni-hybrids.csv"), na.strings = c("")) %>%
  select(
    id, taxon_scientific_name_s_lower, authors_t, rank_s_alphanum,
    # basionym_s_lower, lookup_basionym_id,
    hybrid_parents_s_lower, lookup_hybrid_parent_id
    # , distribution_s_lower
  ) %>%
  mutate(verbatim = paste(taxon_scientific_name_s_lower, authors_t), .keep = "unused") %>%
  mutate(verbatim = str_trim(str_squish(verbatim))) %>%
  # remove taxa with ranks higher than species
  filter(!(rank_s_alphanum %in% c(
    "gen.", "nothosect.", "nothoser.", "nothosubgen.",
    "nothosubsect.", "sect.", "subsect."
  ))) %>%
  # remove parentheses
  mutate(hybrid_parents_s_lower = gsub("^\\(|\\)$", "", hybrid_parents_s_lower)) %>%
  mutate(hybrid_parents_s_lower = str_trim(str_squish(hybrid_parents_s_lower)))

hybrids_ipni <- bind_cols(
  hybrids_ipni %>%
    mutate(parent = str_split(hybrid_parents_s_lower, " × (?=[A-Z])")) %>%
    unnest(cols = "parent"),
  hybrids_ipni %>%
    mutate(parent_taxonID = str_split(lookup_hybrid_parent_id, ",")) %>%
    unnest("parent_taxonID") %>%
    select(parent_taxonID)
) %>%
  distinct() %>%
  select(-rank_s_alphanum, -lookup_hybrid_parent_id, -hybrid_parents_s_lower) %>%
  mutate(parent_taxonID = paste0("urn:lsid:ipni.org:names:", parent_taxonID)) %>%
  mutate(parent_taxonID = ifelse(is.na(parent), NA, parent_taxonID)) %>%
  mutate(edited = verbatim) %>%
  write_csv2(here("data/hybrids_ipni.csv"))
