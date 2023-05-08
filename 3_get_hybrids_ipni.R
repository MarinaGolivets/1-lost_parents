# Marina Golivets
# BiCIKL Hackathon (Sept 20-24, 2021, Meise, Belgium)
# Extract hybrid and parent taxon names from IPNI


# get hybrid names from IPNI -----------------------------------------------------------------------
hybrids_ipni <- data.table::fread(here("input_datasets/hybrids_ipni.csv"), na.strings = c("")) %>%
  select(
    id, taxon_scientific_name_s_lower, authors_t, rank_s_alphanum,
    hybrid_parents_s_lower, lookup_hybrid_parent_id
  ) %>%
  mutate(
    verbatim = paste(taxon_scientific_name_s_lower, authors_t) %>%
      str_squish() %>%
      str_trim(),
    type =  "name",
    .keep = "unused"
  ) %>%
  # remove taxa with ranks higher than species
  filter(!(rank_s_alphanum %in% c(
    "gen.", "nothosect.", "nothoser.", "nothosubgen.",
    "nothosubsect.", "sect.", "subsect."
  ))) %>%
  # remove parentheses
  mutate(
    hybrid_parents_s_lower = str_remove_all(hybrid_parents_s_lower, pattern = "^\\(|\\)$") %>%
      str_squish() %>%
      str_trim()
  )

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
  mutate(
    parent_taxonID = paste0("urn:lsid:ipni.org:names:", parent_taxonID) %>%
      if_else(is.na(parent), NA_character_, .),
    edited = str_replace_all(verbatim, c(" × " = " ×", "^× " = "×")),
    source = "ipni"
  ) %>%
  write_csv2(here("output_datasets/hybrids_ipni.csv"))
