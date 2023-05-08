# Marina Golivets
# BiCIKL Hackathon (Sept 20-24, 2021, Meise, Belgium)
# Extract hybrid and parent taxon names from Ecological Flora of the British Isles (EcoFlora)


# get hybrid names from EcoFlora -------------------------------------------------------------------
# Note: Re-uses the EcoFlora data set was prepared by M. Golivets within a different project
hybrids_ecoflora <- read_csv2(here("input_datasets/taxa_ecoflora.csv")) %>%
  filter(
    is_hybrid == TRUE, # select only hybrids
    str_detect(taxon_name, " '", negate = TRUE) # remove cultivars
  ) %>%
  transmute(
    verbatim = taxon_name,
    type = if_else(str_detect(verbatim, "^[A-Za-z]+ x |^X "), "name", "formula"),
    edited = str_replace_all(verbatim, c(" x " = " × ", "^X " = "×")),
    formula = coalesce(hybrid_formula, taxon_name) %>%
      str_replace_all(c("\\." = ". ", "x jackii" = "× jackii")) %>%
      if_else(str_detect(., "^[A-Za-z]+ x |^X "), NA_character_, .),
    parent = str_split(formula, pattern = " x ")
  ) %>%
  unnest(cols = "parent") %>% # parse hybrid parents from hybrid formulas
  mutate(across(everything(), ~ str_squish(.) %>% str_trim())) %>%
  mutate(
    genus = str_extract(edited, "\\w*"),
    parent = str_replace(parent, "^[A-Z]\\.", genus) %>%
      if_else(str_detect(., "^[a-z]"), paste(genus, .), .) %>%
      if_else(str_detect(., " ", negate = TRUE), NA_character_, .),
    edited = if_else(type == "name", str_replace(edited, " × ", " ×"), edited)
  ) %>%
  select(-genus, -formula) %>%
  mutate(source = "ecoflora") %>%
  distinct() %>%
  glimpse()

w_parent <- hybrids_ecoflora %>%
  filter(!is.na(parent))
wo_parent <- hybrids_ecoflora %>%
  filter(is.na(parent)) %>%
  filter(
    !(edited %in% w_parent$edited),
    edited != "Mentha × villosa Huds."
  )

hybrids_ecoflora <- bind_rows(w_parent, wo_parent) %>%
  write_csv2(here("output_datasets/hybrids_ecoflora.csv"))
