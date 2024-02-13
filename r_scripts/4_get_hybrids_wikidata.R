# Marina Golivets
# BiCIKL Hackathon (Sept 20-24, 2021, Meise, Belgium)
# Extract hybrid and parent taxon names from Wikidata


# load taxon names from Wikidata -------------------------------------------------------------------
hybr_wiki <- readxl::read_excel(
  here("input/Summary_Wikidata.xlsx"),
  sheet = 2, skip = 1
) %>%
  janitor::clean_names() %>%
  transmute(
    verbatim = item_label_2,
    parent = item_label_4
  )

# match to GBIF backbone ---------------------------------------------------------------------------
matched <- gbif.cmpfn(
  taxonName = unique(hybr_wiki$verbatim),
  taxonID = seq_along(unique(hybr_wiki$verbatim))
) %>%
  select(taxonName, scientificName, matchType, confidence, usageKey) %>%
  rename(verbatim = taxonName)

parent_matched <- gbif.cmpfn(
  taxonName = unique(hybr_wiki$parent),
  taxonID = seq_along(unique(hybr_wiki$parent))
) %>%
  select(taxonName, scientificName, usageKey, matchType, confidence) %>%
  rename_with(~ paste0("parent_", .), .cols = everything()) %>%
  mutate(parent = parent_taxonName)

hybr_wiki %<>%
  left_join(matched) %>%
  left_join(parent_matched) %>%
  filter(!is.na(scientificName)) %>%
  mutate(parent = if_else(is.na(parent_taxonName), NA_character_, parent)) %>%
  select(verbatim, parent)

w_parents <- hybr_wiki %>%
  filter(!is.na(parent))
wo_parents <- hybr_wiki %>%
  filter(is.na(parent), !(verbatim %in% w_parents$verbatim))

hybr_wiki <- bind_rows(w_parents, wo_parents) %>%
  mutate(
    type = if_else(str_detect(verbatim, "^[A-Za-z]+ × |^× "), "name", "formula"),
    edited = if_else(
      type == "name", str_replace_all(verbatim, c(" × " = " ×", "^× " = "×")), verbatim
    ),
    type = "name",
    source = "wikidata"
  ) %>%
  distinct() %>%
  write_csv2(here("output/hybrids_wikidata.csv"))
