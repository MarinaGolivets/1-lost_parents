# Marina Golivets
# Extract hybrid parentage from:
# 1. IPNI
# 2. World Flora Online (WFO)


# LOAD WCVP DATA -----------------------------------------------------------------------------------

wcvp_v11 <- fread(here("wcvp_v11/wcvp_names.csv"), quote = "", encoding = "UTF-8") %>%
  mutate(
    name_to_match = str_remove_all(taxon_name, "× |-|\\+") %>%
      str_c(taxon_authors, sep = " ") %>%
      str_squish() %>%
      str_trim() %>%
      stri_trans_general("latin-ascii"),
    across(c(taxon_rank, taxon_status), ~ str_to_sentence(.)),
    across(where(is.character), ~ na_if(., "")),
    source = "WCVP",
    ver = "v11"
  ) %>%
  glimpse()


# ADD HYBRID PARENTAGE INFO FROM WCVP --------------------------------------------------------------

all_h2 <- read_rds(here("outputs/hybrids-wfo-and-powo_0.rds")) %>%
  mutate(acc_local_id = coalesce(acc_local_id, local_id)) %>%
  left_join(select(wcvp_v11, powo_id, hybrid_formula)) %>%
  group_by(ipni_id) %>%
  fill(hybrid_formula, .direction = "downup") %>%
  ungroup() %>%
  mutate(hybrid_formula = str_remove_all(hybrid_formula, "^\\(|\\)$"))


# ADD HYBRID PARENTAGE INFO FROM IPNI --------------------------------------------------------------

# set up parallel computing
cl <- makeCluster(detectCores() - 1)
doParallel::registerDoParallel(cl)

# pull IPNI IDs
ipni_ids <- unique(all_h2$ipni_id)
ipni_ids <- ipni_ids[!is.na(ipni_ids)]

# search IPNI
start_time <- Sys.time()
q_id_ipni_l <- foreach(
  i = seq_along(ipni_ids),
  .packages = c("kewr"),
  .combine = "c",
  .errorhandling = "pass"
) %dopar% {
  out <- lookup_ipni(ipni_ids[i], type = "taxon")
  return(list(out))
}
end_time <- Sys.time()
end_time - start_time
write_rds(q_id_ipni_l, here("q2_id_ipni_l.rds"))

# unlist IPNI query results
ipni_parents <- q_id_ipni_l %>%
  .[sapply(., typeof) == "list"] %>%
  map(., ~ .[c("id", "hybridParents")] %>% compact()) %>%
  data.table::rbindlist(., fill = TRUE) %>%
  distinct() %>%
  rename(ipni_id = id) %>%
  mutate(parent_seq = seq_along(ipni_id), .by = ipni_id) %>%
  unnest_longer(hybridParents) %>%
  filter(hybridParents_id %in% c("id", "wfoId", "name", "authors")) %>%
  mutate(hybridParents = unlist(hybridParents)) %>%
  pivot_wider(names_from = hybridParents_id, values_from = hybridParents) %>%
  rename(
    parent_name = name, parent_authors = authors,
    parent_ipni_id = id, parent_wfo_id = wfoId
  ) 


ipni_formula <- q_id_ipni_l %>%
  .[sapply(., typeof) == "list"] %>%
  map(., ~ .[c("id", "originalHybridParentage", "remarks", "originalRemarks")] %>%
    compact()) %>%
  data.table::rbindlist(., fill = TRUE) %>%
  distinct() %>%
  rename(ipni_id = id) %>%
  mutate(
    hybrid_formula = str_remove_all(originalHybridParentage, "^\\(|\\)\\.$") %>%
      coalesce(str_extract(originalRemarks, "\\(.* [xX] .*\\)") %>%
        str_remove_all("^\\(.*\\(|^\\(|\\)$")) %>%
      coalesce(str_extract(remarks, ".* [xX] .*")) %>%
      str_remove(
        "^.*[Pp]arents: |^.*parents are: |^.*may be: |^.*parents \\?: |^.*parents \\(\\?\\):"),
    artificial_hybrid = str_detect(remarks, "[cC]ult|[hH]ort|garden|artef") %>%
      na_if(FALSE) %>%
      coalesce(str_detect(originalRemarks, "[cC]ult|[hH]ort|garden|artef")) %>%
      na_if(FALSE)
  )

ipni_formula <- ipni_parents %>%
  group_by(ipni_id) %>%
  summarise(hybrid_formula = str_c(parent_name, collapse = " × ")) %>%
  filter(str_detect(hybrid_formula, " × ")) %>%
  full_join(ipni_formula, by = c("ipni_id" = "ipni_id")) %>%
  mutate(hybrid_formula = coalesce(hybrid_formula.y, hybrid_formula.x), .keep = "unused") 


all_h3 <- all_h2 %>%
  select(ipni_id, hybrid_formula) %>%
  filter(!is.na(ipni_id)) %>%
  full_join(
    select(ipni_formula, ipni_id, hybrid_formula, artificial_hybrid) %>%
      filter(!is.na(hybrid_formula)),
      relationship = "many-to-many", na_matches = "never"
  ) %>%
  right_join(select(all_h2, -hybrid_formula), relationship = "many-to-many") %>%
  mutate(hybrid_formula = str_replace_all(hybrid_formula, " x | X ", " × ") %>%
           str_squish() %>%
           str_trim()) %>%
  mutate(
    genus = str_extract(name_to_match, "\\w*"),
    hybrid_formula = str_replace(hybrid_formula, "^[A-Z]\\.", genus) %>%
      str_replace("^\\? [A-Z]\\.", str_c("? ", genus)) %>%
      str_replace_all(" × [A-Z]\\.", str_c(" × ", genus)),
    artificial_hybrid = if_else(taxon_status == "Artificial hybrid", TRUE, artificial_hybrid)
  ) %>%
  select(-genus) %>%
  distinct() %>%
  group_by(name_to_match) %>%
  fill(c(hybrid_formula, artificial_hybrid), .direction = "downup") %>%
  group_by(ipni_id) %>%
  fill(c(hybrid_formula, artificial_hybrid), .direction = "downup") %>%
  ungroup()


# keep only IPNI top copy 
top_copy <- q_id_ipni_l %>%
  .[sapply(., typeof) == "list"] %>%
  map(., ~ .[c("id", "topCopy", "suppressed")] %>% 
        compact()) %>%
  data.table::rbindlist(., fill = TRUE) %>%
  rename(ipni_id = id) %>%
  clean_names() %>%
  distinct()

all_h4 <- all_h3 %>%
  left_join(top_copy) %>%
  filter(top_copy == TRUE | is.na(top_copy)) %>%
  filter(n() == 1 | n() > 1 & suppressed == FALSE, .by = local_id) %>%
  distinct()
 
all_h4 %<>%
  select(local_id, ipni_id) %>%
  distinct() %>%
  group_by(local_id) %>%
  summarise(ipni_id = str_c(ipni_id, collapse = "|")) %>%
  right_join(select(all_h4, -ipni_id, -top_copy, -suppressed)) %>%
  distinct()

# ADD HYBRID PARENTAGE INFO FROM VARIOUS SOURCES ---------------------------------------------------

gbif <- read_csv2(here("outputs/hybrids_gbif_checklists.csv")) %>%
  filter(type == "name") %>%
  transmute(full_name = edited, parent) %>%
  distinct() %>%
  group_by(full_name) %>%
  summarise(hybrid_formula = str_c(parent, collapse = " × "), .groups = "drop") %>%
  mutate(
    name_to_match = str_remove_all(full_name, "×|-") %>%
      stri_trans_general("latin-ascii"),
    source = "GBIF"
  )

gbif_occ <- read_csv2(here("outputs/hybrids_gbif_occurrences.csv")) %>%
  filter(type == "name", !is.na(parent)) %>%
  transmute(full_name = edited, parent) %>%
  distinct() %>%
  group_by(full_name) %>%
  summarise(hybrid_formula = str_c(parent, collapse = " × "), .groups = "drop") %>%
  mutate(
    name_to_match = str_remove_all(full_name, "×|-") %>%
      stri_trans_general("latin-ascii"),
    source = "GBIF"
  )

ecoflora <- read_csv2(here("outputs/hybrid-formulas_ecoflora.csv")) %>%
  filter(type == "name") %>%
  transmute(
    full_name = edited,
    name_to_match = str_remove_all(edited, "×|-") %>%
      str_squish() %>%
      str_trim() %>%
      stri_trans_general("latin-ascii"),
    hybrid_formula = formula
  ) %>%
  mutate(source = "EcoFlora")

floras <- read_csv2(here("outputs/hybrids_floras.csv")) %>%
  filter(type == "name") %>%
  transmute(full_name = edited, source)

wikidata <- read_csv2(here("outputs/hybrids_wikidata.csv")) %>%
  filter(type == "name") %>%
  transmute(full_name = edited, parent) %>%
  distinct() %>%
  group_by(full_name) %>%
  summarise(hybrid_formula = str_c(parent, collapse = " × "), .groups = "drop") %>%
  mutate(
    name_to_match = str_remove_all(full_name, "×|-") %>%
      stri_trans_general("latin-ascii"),
    source = "Wikidata"
  )

all_h5 <- all_h4 %>% 
  mutate(full_name = if_else(
    is.na(taxon_authors), taxon_name, str_c(taxon_name, " ", taxon_authors))
    ) %>%
  bind_rows(gbif, gbif_occ, ecoflora, floras, wikidata) %>%
  distinct() %>%
  group_by(name_to_match) %>%
  fill(c(hybrid_formula, artificial_hybrid), .direction = "downup") %>%
  select(-taxon_name, -taxon_authors, -scientific_name_id) %>%
  relocate(
    local_id, acc_local_id, ipni_id, powo_id, wfo_id,
    full_name, name_to_match, taxon_rank,	taxon_status,	nomenclatural_status,
    artificial_hybrid, hybrid_formula
    ) %>%
  arrange(name_to_match)

writexl::write_xlsx(all_h5, here("all-hybrids_in-work.xlsx"))

all_h5 <- readxl::read_xlsx(here("all-hybrids_in-work.xlsx"))
# pull IPNI IDs
powo_ids <- unique(all_h5$powo_id)
powo_ids <- powo_ids[!is.na(powo_ids)]

# search IPNI
start_time <- Sys.time()
q_id_powo_l <- foreach(
  i = seq_along(powo_ids),
  .packages = c("kewr"),
  .combine = "c",
  .errorhandling = "pass"
) %dopar% {
  out <- lookup_powo(powo_ids[i])
  return(list(out))
}
end_time <- Sys.time()
end_time - start_time

write_rds(q_id_powo_l, here("q_id_powo_l.rds"))
