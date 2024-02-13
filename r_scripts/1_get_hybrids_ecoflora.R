# Marina Golivets
# BiCIKL Hackathon (Sept 20-24, 2021, Meise, Belgium)
# Extract hybrid and parent taxon names from Ecological Flora of the British Isles (EcoFlora)

# Note: Here we reuse the EcoFlora data set was prepared by M. Golivets within a different project

# get hybrid names from EcoFlora -------------------------------------------------------------------
hybr_ecoflora <- read_csv2(here("input/taxa_ecoflora.csv")) %>%
  filter(
    is_hybrid == TRUE, # select only hybrids
    str_detect(taxon_name, " '| hort\\.|Astilbe x arendsii", negate = TRUE) # remove cultivars
  ) %>%
  transmute(
    verbatim = taxon_name,
    type = if_else(str_detect(verbatim, "^[A-Za-z]+ x |^X "), "name", "formula"),
    edited = str_replace_all(
      verbatim, c
      (" x " = " × ", "^X " = "×", "Qercus" = "Quercus", " maueri " = " maureri ")
    ),
    formula = coalesce(hybrid_formula, taxon_name) %>%
      str_replace_all(
        c("\\." = ". ", "x jackii" = "× jackii", "niedzwetskyana" = "niedzwetzkyana")) %>%
      if_else(str_detect(., "^[A-Za-z]+ x |^X "), NA_character_, .),
    parent = str_split(formula, pattern = " x ")
  ) %>%
  unnest(cols = "parent") %>% # parse out hybrid parents from hybrid formulas
  mutate(across(everything(), ~ str_squish(.) %>% str_trim())) %>%
  mutate(
    genus = str_extract(edited, "\\w*"),
    parent = str_replace(parent, "^[A-Z]\\.", genus) %>%
      if_else(str_detect(., "^[a-z]"), paste(genus, .), .) %>%
      if_else(str_detect(., " ", negate = TRUE), NA_character_, .),
    edited2 = if_else(type == "name", str_replace(edited, " × ", " ×"), edited),
    edited3 = str_extract(edited, ".* × \\w*|^× \\w* \\w*") %>%
      if_else(type == "formula", NA_character_, .)
  ) %>%
  filter(str_detect(edited, "× .* ×", negate = TRUE)) %>%
  select(-genus, -formula) %>%
  mutate(source = "ecoflora") %>%
  distinct() %>%
  glimpse()

# remove duplicates
w_parent <- hybr_ecoflora %>%
  filter(!is.na(parent))
wo_parent <- hybr_ecoflora %>%
  filter(is.na(parent)) %>%
  filter(
    !(edited %in% w_parent$edited),
    edited != "Mentha × villosa Huds."
  )
hybr_ecoflora <- bind_rows(w_parent, wo_parent)
rm(w_parent, wo_parent)


# match names to GBIF ------------------------------------------------------------------------------
gbif_h <- gbif.cmpfn(
  unique(hybr_ecoflora$edited2),
  seq_along(unique(hybr_ecoflora$edited2))
) %>%
  filter(str_detect(scientificName, "×| x ") | confidence > 95 | confidence == 0)
setdiff(unique(hybr_ecoflora$edited2), gbif_h$name)

hybr_ecoflora %<>% 
  left_join(
    gbif_h %>% 
    transmute(
      name, usageKey,
      acceptedUsageKey = coalesce(acceptedUsageKey, usageKey)),
    by = c("edited2" = "name")) 

gbif_p <- gbif.cmpfn(
  unique(hybr_ecoflora$parent),
  seq_along(unique(hybr_ecoflora$parent))
) %>%
  filter(str_detect(scientificName, "×| x ") | confidence > 95 | confidence == 0) %>%
  group_by(name) %>%
  filter(n() == 1) %>% 
  ungroup()
setdiff(unique(hybr_ecoflora$parent), gbif_p$name)

hybr_ecoflora %<>% 
  group_by(edited) %>%
  filter(n() == 2) %>%
  mutate(parent_type = c("parent1", "parent2")) %>%
  ungroup() %>%
  pivot_wider(values_from = parent, names_from = parent_type)

hybr_ecoflora %<>% 
  left_join(
    gbif_p %>% 
      transmute(
        parent1 = name, 
        parent1_usageKey = usageKey,
        parent1_acceptedUsageKey = coalesce(acceptedUsageKey, usageKey))
    ) %>%
  left_join(
    gbif_p %>% 
      transmute(
        parent2 = name, 
        parent2_usageKey = usageKey,
        parent2_acceptedUsageKey = coalesce(acceptedUsageKey, usageKey))
  ) 

# save data
write_csv2(hybr_ecoflora, here("output/hybrids_ecoflora.csv"))
