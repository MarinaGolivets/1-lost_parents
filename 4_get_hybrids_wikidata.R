# Marina Golivets
# BiCIKL Hackathon (Sept 20-24, 2021, Meise, Belgium)
# Extract hybrid and parent taxon names from Wikidata


# load taxon names from Wikidata -------------------------------------------------------------------
hybrids_wiki <- readxl::read_excel(
  here("input_datasets/Summary_Wikidata.xlsx"),
  sheet = 2, skip = 1
) %>%
  janitor::clean_names() %>%
  transmute(verbatim = item_label_2, parent = item_label_4)

# match to GBIF backbone ---------------------------------------------------------------------------
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
  mutate(parent = if_else(is.na(parent_taxonName), NA_character_, parent)) %>%
  select(verbatim, parent)

w_parents <- hybrids_wiki %>%
  filter(!is.na(parent))
wo_parents <- hybrids_wiki %>%
  filter(is.na(parent), !(verbatim %in% w_parents$verbatim))

hybrids_wiki <- bind_rows(w_parents, wo_parents) %>%
  mutate(
    type = if_else(str_detect(verbatim, "^[A-Za-z]+ × |^× "), "name", "formula"),
    edited = if_else(
      type == "name", str_replace_all(verbatim, c(" × " = " ×", "^× " = "×")), verbatim),
    type = "name",
    source = "wikidata"
  ) %>%
  distinct() %>%
  write_csv2(here("output_datasets/hybrids_wikidata.csv"))