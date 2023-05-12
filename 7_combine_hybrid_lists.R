# BiCIKL Hackathon (Sept 20-24, 2021, Meise)
# Combine all hybrid lists and match taxon names to the GBIF backbone


# load data sets -----------------------------------------------------------------------------------
gbif_occ <- read_csv2(here("output_datasets/hybrids_gbif_occurrences.csv"))
gbif_chl <- read_csv2(here("output_datasets/hybrids_gbif_checklists.csv"))
wiki <- read_csv2(here("output_datasets/hybrids_wikidata.csv"))
ecoflora <- read_csv2(here("output_datasets/hybrids_ecoflora.csv"))
ipni <- read_csv2(here("output_datasets/hybrids_ipni.csv"))
floras <- read_csv2(here("output_datasets/hybrids_floras.csv"))
vascan <- tibble(
  verbatim = read_lines(here("input_datasets/vascan_hybrid_names.csv"))
  ) %>%
  mutate(
    edited = str_replace_all(verbatim, "×", " × ") %>%
      str_squish() %>%
      str_trim()
)

hybrids0 <- bind_rows(gbif_occ, gbif_chl, wiki, ecoflora, ipni, floras, vascan) %>%
  select(-id, -parent_taxonID) %>%
  distinct()

# match to GBIF ------------------------------------------------------------------------------------
matched0 <- gbif.cmpfn(
  taxonName = unique(hybrids0$edited),
  taxonID = seq_along(unique(hybrids0$edited))
) 

matched <- matched0 %>%
  transmute(
    edited = taxonName, 
    usageKey = coalesce(usageKey, acceptedUsageKey),
    confidence, matchType, scientificName
    ) %>%
  filter(
    !(str_detect(scientificName, " subsp\\.") &
    str_detect(edited, "subsp\\.|var\\.") &
      matchType == "FUZZY"),
    !(str_detect(scientificName, " var\\.|×") &
    str_detect(edited, "var\\.|subsp\\.") &
      matchType == "FUZZY")
    ) %>%
  distinct()
  
# match parent taxa
parent_matched0 <- gbif.cmpfn(
  taxonName = unique(hybrids0$parent),
  taxonID = seq_along(unique(hybrids0$parent))
)

parent_matched <- parent_matched0 %>%
  transmute(
    taxonName, 
    usageKey = coalesce(usageKey, acceptedUsageKey),
    confidence, matchType, scientificName
  ) %>%
  rename_with(~ paste0("parent_", .), .cols = everything()) %>%
  rename(parent = parent_taxonName) %>%
  distinct()

hybrids_matched <- hybrids0 %>%
  left_join(matched) %>%
  left_join(parent_matched) %>%
  filter(!is.na(usageKey)) %>%
  select(-verbatim, -parent) %>%
  distinct()

hybrids <- hybrids0 %>%
  left_join(matched) %>%
  left_join(parent_matched) %>%
  filter(is.na(usageKey)) %>%
  bind_rows(hybrids_matched) %>%
  select(-verbatim) %>%
  distinct()

w_parents <- hybrids %>%
  filter(!is.na(parent_scientificName))
wo_parents <- hybrids %>%
  filter(is.na(parent_scientificName)) %>%
  filter(!(scientificName %in% w_parents$scientificName))
hybrids <- bind_rows(w_parents, wo_parents)

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
  mutate(parent = coalesce(parent_scientificName, parent),
         hybrid = coalesce(scientificName, edited)) %>%
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
  transmute(hybrid = coalesce(scientificName, edited), 
            parent, parent_name) %>%
  pivot_wider(values_from = parent, names_from = parent_name) 

# exclude hybrids with more than one name
hybrids2 %<>% group_by(parent1, parent2) %>% 
  filter(n() == 1) %>%
  ungroup()

hybrids2 %<>% 
  left_join(
  hybrids_matched %>%
    select(usageKey, scientificName),
  by = c("hybrid" = "scientificName")) %>%
  left_join(
    parent_matched %>%
      select(parent_usageKey, parent_scientificName),
    by = c("parent1" = "parent_scientificName")) %>%
  rename(parent_usageKey1 = parent_usageKey) %>%
  left_join(
    parent_matched %>%
      select(parent_usageKey, parent_scientificName),
    by = c("parent2" = "parent_scientificName")) %>%
  rename(parent_usageKey2 = parent_usageKey) %>%
  distinct() %>%
  filter(str_detect(hybrid, " × ", negate = TRUE))

write_csv2(hybrids2, here("output_datasets/hybrids_with_two_parents.csv"))
