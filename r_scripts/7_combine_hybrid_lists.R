# BiCIKL Hackathon (Sept 20-24, 2021, Meise)
# Combine all hybrid lists and match taxon names to the GBIF backbone


# load data sets -----------------------------------------------------------------------------------
ipni <- read_csv2(here("output/hybrids_ipni.csv"))
gbif_occ <- read_csv2(here("output/hybrids_gbif_occurrences.csv"))
gbif_chl <- read_csv2(here("output/hybrids_gbif_checklists.csv"))
wiki <- read_csv2(here("output/hybrids_wikidata.csv"))
ecoflora <- read_csv2(here("output/hybrids_ecoflora.csv"))
floras <- read_csv2(here("output/hybrids_floras.csv"))
vascan <- tibble(
  verbatim = read_lines(here("input/vascan_hybrid_names.csv")),
  source = "vascan"
)

# combine data sets
hybrids0 <- bind_rows(gbif_occ, gbif_chl, wiki, ecoflora, floras, vascan, ipni) %>%
  select(-id, -parent_taxonID, -verbatim) %>%
  mutate(edited = str_replace(edited, " × [A-Z]\\. ", " × ")) %>%
  distinct()

# match hybrid names to GBIF -----------------------------------------------------------------------
matched0 <- gbif.cmpfn(
  taxonName = unique(hybrids0$edited),
  taxonID = seq_along(unique(hybrids0$edited))
)

matched <- matched0 %>%
  transmute(
    edited = taxonName,
    usageKey = coalesce(acceptedUsageKey, usageKey),
    matchType, confidence, scientificName
  ) %>%
  filter(
    confidence >= 75,
    !(str_detect(edited, " subsp\\.|var\\.", negate = TRUE) &
      str_detect(scientificName, "subsp\\.|var\\.") &
      matchType == "FUZZY"),
    !(str_detect(scientificName, " var\\.|×|subsp\\.") &
      str_detect(edited, "var\\.|subsp\\.") &
      matchType == "FUZZY")
  ) %>%
  distinct() %>%
  group_by(edited) %>%
  filter(n() == 1) %>%
  ungroup()
setdiff(matched0$taxonName, matched1$edited)

hybrids0 %<>%
  left_join(matched)
hybrids1 <- hybrids0 %>%
  mutate(
    scientificName = coalesce(scientificName, edited), 
    .keep = "unused") %>%
  select(-matchType, -confidence) %>%
  distinct() %>%
  group_by(scientificName, usageKey, parent, type) %>%
  summarise(source = str_c(source, collapse = "|"), .groups = "drop")


w_parents <- hybrids1 %>%
  filter(!is.na(parent))
wo_parents <- hybrids1 %>%
  filter(is.na(parent)) %>%
  filter(!(usageKey %in% w_parents$usageKey))
hybrids1 <- bind_rows(w_parents, wo_parents) %>%
  group_by(edited) %>%
  filter(n() < 3) # exclude hybrids w/ > 2 parents

# ipni_ <- hybrids1 %>%
#   filter(!is.na(parent) & source == "ipni" & !is.na(usageKey))
# not_ipni_ <- hybrids1 %>%
#   filter(!(usageKey %in% ipni_$usageKey))
# hybrids2 <- bind_rows(ipni_, not_ipni_)


# match parent taxa to GBIF ------------------------------------------------------------------------
parent_matched0 <- gbif.cmpfn(
  taxonName = unique(hybrids1$parent),
  taxonID = seq_along(unique(hybrids1$parent))
)

parent_matched <- parent_matched0 %>%
  filter(confidence >= 75) %>%
  transmute(
    taxonName,
    usageKey = coalesce(acceptedUsageKey, usageKey),
    confidence, matchType, scientificName
  ) %>%
  rename_with(~ paste0("parent_", .), .cols = everything()) %>%
  rename(parent = parent_taxonName) %>%
  group_by(parent) %>%
  filter(n() == 1) %>%
  ungroup() %>%
  distinct()

hybrids1 %<>%
  left_join(parent_matched)
hybrids2 <- hybrids1 %>%
  transmute(
    scientificName = coalesce(scientificName, edited),
    parent_scientificName = coalesce(parent_scientificName, parent),
    source
  ) %>%
  distinct() %>%
  group_by(scientificName, parent_scientificName) %>%
  summarise(source = str_c(source, collapse = "|"), .groups = "drop")

  
  group_by(scientificName, source) %>%
  arrange(parent_scientificName) %>%
  filter(n() == 2) %>%
  mutate(parent_name = c("parent1", "parent2")) %>%
  ungroup() %>%
  pivot_wider(values_from = parent_scientificName, names_from = parent_name)

names_ <- hybrids2 %>%
  group_by(parent1, parent2) %>%
  mutate(group_id = cur_group_id()) %>%
  filter(type == "name")
formulas_ <- hybrids2 %>%
  group_by(parent1, parent2) %>%
  mutate(group_id = cur_group_id()) %>%
  filter(!(group_id %in% names_$group_id))
hybrids3 <- bind_rows(names_, formulas_)

hybrids6 <- hybrids5 %>%
  left_join(
    hybrids3 %>%
      select(scientificName, usageKey) %>%
      distinct()
  ) %>%
  left_join(
    hybrids3 %>%
      select(parent_scientificName, parent_usageKey) %>%
      rename(parent1_usageKey = parent_usageKey) %>%
      distinct(),
    by = c("parent1" = "parent_scientificName")
  ) %>%
  left_join(hybrids3 %>%
    select(parent_scientificName, parent_usageKey) %>%
    rename(parent2_usageKey = parent_usageKey) %>%
    distinct(),
  by = c("parent2" = "parent_scientificName")
  ) %>%
  group_by(parent1_usageKey, parent2_usageKey) %>%
  mutate(group_id = cur_group_id())

names_ <- hybrids6 %>%
  filter(type == "name")
formulas_ <- hybrids6 %>%
  filter(!(group_id %in% names_$group_id))
hybrids7 <- bind_rows(names_, formulas_)

groups_resolved <- hybrids6 %>%
  filter(n() > 1 & !is.na(usageKey))
hybrids8 <- hybrids7 %>%
  filter(!(group_id %in% groups_resolved$group_id & is.na(usageKey)))


matched_accepted <- pbapply::pblapply(unique(hybrids8$usageKey), get_acc_name_gbif.fn, cl = cl)
gbif_acc <- gbif_acc0 %>% bind_rows()

hybrids8 <- hybrids8 %>%
  select(-source) %>%
  left_join(
    hybrids8 %>%
      group_by(group_id) %>%
      summarise(source = str_c(source, collapse = "|"))
  ) %>%
  distinct()

n_distinct(hybrids$usageKey) # 20999
n_distinct(hybrids %>% filter(is.na(usageKey)) %>% select(edited)) # 25182

# matched hybrids with matched parents
hybrids %>%
  filter(!is.na(usageKey) & !is.na(parent_usageKey)) %>%
  select(usageKey, parent_usageKey) %>%
  distinct() %>%
  pull(usageKey) %>%
  n_distinct() # 6496

# non-matched hybrids with parents
hybrids %>%
  filter(is.na(usageKey) & !is.na(parent)) %>%
  select(edited, parent) %>%
  distinct() %>%
  pull(edited) %>%
  n_distinct() # 25182

hybrids1 <- hybrids %>%
  mutate(
    parent = coalesce(parent_scientificName, parent),
    hybrid = coalesce(scientificName, edited)
  ) %>%
  group_by(usageKey, parent) %>%
  filter(confidence == max(confidence) | is.na(confidence)) %>%
  filter(n() == 1) %>%
  group_by(hybrid) %>%
  filter(n() > 1)

# select hybrids with two parents
hybrids2 <- hybrids1 %>%
  arrange(hybrid, parent) %>%
  filter(n() == 2) %>%
  mutate(parent_name = c("parent1", "parent2")) %>%
  ungroup() %>%
  transmute(
    hybrid = coalesce(scientificName, edited),
    parent, parent_name
  ) %>%
  pivot_wider(values_from = parent, names_from = parent_name)

# exclude hybrids with more than one name
hybrids2 %<>% group_by(parent1, parent2) %>%
  filter(n() == 1) %>%
  ungroup()

hybrids2 %<>%
  left_join(
    hybrids_matched %>%
      select(usageKey, scientificName),
    by = c("hybrid" = "scientificName")
  ) %>%
  left_join(
    parent_matched %>%
      select(parent_usageKey, parent_scientificName),
    by = c("parent1" = "parent_scientificName")
  ) %>%
  rename(parent_usageKey1 = parent_usageKey) %>%
  left_join(
    parent_matched %>%
      select(parent_usageKey, parent_scientificName),
    by = c("parent2" = "parent_scientificName")
  ) %>%
  rename(parent_usageKey2 = parent_usageKey) %>%
  distinct() %>%
  filter(str_detect(hybrid, " × ", negate = TRUE))

write_csv2(hybrids2, here("output/hybrids_with_two_parents.csv"))
