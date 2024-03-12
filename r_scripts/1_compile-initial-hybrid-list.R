# Marina Golivets
# Extract all hybrid names (accepted and synonyms) from:
# 1. World Flora Online (WFO)
# 2. World Checklist of Vascular Plants (WCVP)


# **************************************************************************************************
# LOAD WFO DATA ------------------------------------------------------------------------------------

# NOTE: I use several latest versions of WFO because some hybrids are excluded
# and some are unidentifiable in the newer versions

wfo_v1.1 <- fread(
  here("WFOTaxonomicBackbone_v.1.1_20230107/classification.txt"),
  select = c(
    "taxonID", "acceptedNameUsageID", "scientificNameID", "scientificName",
    "scientificNameAuthorship", "taxonRank", "taxonomicStatus", 
    "nomenclaturalStatus", "source", "taxonRemarks"
  ),
  quote = "", encoding = "UTF-8"
) %>%
  mutate(
    across(where(is.character), ~ str_replace_all(., c('^"|"$' = ""))),
    wfo_id = taxonID,
    name_to_match = str_remove_all(scientificName, "×|-") %>%
      str_c(scientificNameAuthorship, sep = " ") %>%
      str_squish() %>%
      str_trim() %>%
      stri_trans_general("latin-ascii"),
    across(c(taxonRank, taxonomicStatus), ~ str_to_sentence(.)),
    ipni_id = str_remove(scientificNameID, "urn:lsid:ipni.org:names:"),
    across(where(is.character), ~ na_if(., "")),
    original_source = source,
    source = "WFO",
    ver = "v1.1"
  ) %>%
  clean_names() %>%
  glimpse()


wfo_v2.1 <- fread(
  here("WFOTaxonomicBackbone_v.2.1_20230304/classification.csv"),
  select = c(
    "taxonID", "acceptedNameUsageID", "scientificNameID", "scientificName",
    "scientificNameAuthorship", "taxonRank", "taxonomicStatus", 
    "nomenclaturalStatus", "source", "taxonRemarks"
  ),
  quote = "", encoding = "UTF-8"
) %>%
  # clean raw taxon names
  mutate(scientificName = map_chr(
    str_extract_all(scientificName, "[A-Za-z-×öëï. ]+"), ~ str_c(., collapse = "")
  )) %>%
  mutate(
    across(where(is.character), ~ str_replace_all(., c('^"|"$' = ""))),
    wfo_id = taxonID,
    name_to_match = str_remove_all(scientificName, "×|-") %>%
      str_c(scientificNameAuthorship, sep = " ") %>%
      str_squish() %>%
      str_trim() %>%
      stri_trans_general("latin-ascii"),
    across(c(taxonRank, taxonomicStatus), ~ str_to_sentence(.)),
    ipni_id = str_remove(scientificNameID, "urn:lsid:ipni.org:names:") %>%
      if_else(str_detect(., "-"), ., NA_character_),
    across(where(is.character), ~ na_if(., "")),
    original_source = source,
    source = "WFO",
    ver = "v2.1"
  ) %>%
  clean_names() %>%
  glimpse()


wfo_v3 <- fread(
  here("wfo_202312.csv"),
  select = c(
    "taxonID", "acceptedNameUsageID", "scientificNameID", "scientificName",
    "scientificNameAuthorship", "taxonRank", "taxonomicStatus", 
    "nomenclaturalStatus", "source", "taxonRemarks"
  ),
  quote = "", encoding = "UTF-8"
) %>%
  mutate(
    across(where(is.character), ~ str_replace_all(., c('^"|"$' = ""))),
    wfo_id = taxonID,
    name_to_match = str_remove_all(scientificName, "×|-") %>%
      str_c(scientificNameAuthorship, sep = " ") %>%
      str_squish() %>%
      str_trim() %>%
      stri_trans_general("latin-ascii"),
    across(c(taxonRank, taxonomicStatus), ~ str_to_sentence(.)),
    ipni_id = str_remove(scientificNameID, "urn:lsid:ipni.org:names:") %>%
      if_else(str_detect(., "-"), ., NA_character_),
    across(where(is.character), ~ na_if(., "")),
    original_source = source,
    source = "WFO",
    ver = "v2023-12"
  ) %>%
  clean_names() %>%
  glimpse()


# PROVIDE ADDITIONAL IPNI IDS FOR NAMES IN WFO -----------------------------------------------------

# fill in additional IPNI IDs (from Zenodo)
wfo_to_ipni <- fread(here("inputs/ipni_to_wfo_v2.1.csv")) %>%
  bind_rows(fread(here("inputs/ipni_to_wfo_202312.csv"))) %>%
  mutate(ipni_id = str_remove(ipni_id, "urn:lsid:ipni.org:names:")) %>%
  distinct()

wfo_v1.1 <- wfo_to_ipni %>%
  inner_join(select(wfo_v1.1, -ipni_id)) %>%
  full_join(wfo_v1.1) %>%
  filter(n() == 1 | (n() > 1 & !is.na(ipni_id)), .by = taxon_id)

wfo_v2.1 <- wfo_to_ipni %>%
  inner_join(select(wfo_v2.1, -ipni_id)) %>%
  full_join(wfo_v2.1) %>%
  filter(n() == 1 | (n() > 1 & !is.na(ipni_id)), .by = taxon_id)

wfo_v3 <- wfo_to_ipni %>%
  inner_join(select(wfo_v3, -ipni_id)) %>%
  full_join(wfo_v3) %>%
  filter(n() == 1 | (n() > 1 & !is.na(ipni_id)), .by = taxon_id)


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


# SUBSET AN INITIAL LIST OF HYBRIDS FROM WCVP ------------------------------------------------------

# get hybrids from the backbones based on the occurrence of "×" in the taxon names

wcvp_v11_h0 <- wcvp_v11 %>%
  filter(str_detect(taxon_name, "×|\\+"), taxon_rank != "Genus")
# include accepted names
wcvp_v11_h0 <- wcvp_v11 %>%
  filter(plant_name_id %in% c(wcvp_v11_h0$plant_name_id, wcvp_v11_h0$accepted_plant_name_id))


# SUBSET AN INITIAL LIST OF HYBRIDS FROM WFO -------------------------------------------------------

wfo_v1.1_h0 <- wfo_v1.1 %>%
  filter(
    str_detect(scientific_name, "×") | name_to_match %in% wcvp_v11_h0$name_to_match |
      !is.na(ipni_id) & ipni_id %in% wcvp_v11_h0$ipni_id,
    taxon_rank != "Genus"
  )

wfo_v2.1_h0 <- wfo_v2.1 %>%
  filter(
    str_detect(scientific_name, "×") | taxon_id %in% wfo_v1.1_h0$taxon_id |
      name_to_match %in% wcvp_v11_h0$name_to_match | 
      !is.na(ipni_id) & ipni_id %in% wcvp_v11_h0$ipni_id |
      str_detect(taxon_remarks, "[hH]ybr|×"),
    taxon_rank != "Genus"
  )

wfo_v3_h0 <- wfo_v3 %>%
  filter(
    str_detect(scientific_name, "×") |
      taxon_id %in% c(wfo_v1.1_h0$taxon_id, wfo_v2.1_h0$taxon_id) |
      name_to_match %in% wcvp_v11_h0$name_to_match |
      !is.na(ipni_id) & ipni_id %in% wcvp_v11_h0$ipni_id |
      str_detect(taxon_remarks, "[hH]ybr|×"),
    taxon_rank != "Genus"
  )

# (from earlier versions of WFO use only the names that are not in the latest version)
wfo_v1.1_h0 %<>%
  filter(!(taxon_id %in% c(wfo_v2.1_h0$taxon_id, wfo_v3_h0$taxon_id)))
wfo_v2.1_h0 %<>%
  filter(!(taxon_id %in% wfo_v3_h0$taxon_id))

# combine hybrids across WFO versions
wfo_h0 <- bind_rows(
  wfo_v1.1 %>%
    # subset hybrids and their accepted names (accepted names are not hybrids sometimes)
    filter(taxon_id %in% c(wfo_v1.1_h0$taxon_id, wfo_v1.1_h0$accepted_name_usage_id)),
  wfo_v2.1 %>%
    filter(
      taxon_id %in% c(
        wfo_v2.1_h0$taxon_id, wfo_v2.1_h0$accepted_name_usage_id,
        wfo_v1.1_h0$accepted_name_usage_id
      )
    )
) %>%
  filter(n() == 1 | (n() > 1 & ver == "v2.1"), .by = taxon_id) %>%
  bind_rows(
    wfo_v3 %>%
      filter(
        taxon_id %in% c(
          wfo_v3_h0$taxon_id, wfo_v3_h0$accepted_name_usage_id,
          wfo_v2.1_h0$accepted_name_usage_id, wfo_v1.1_h0$accepted_name_usage_id
        )
      )
  ) %>%
  filter(n() == 1 | (n() > 1 & ver == "v2023-12"), .by = taxon_id) %>%
  mutate(
    local_id = taxon_id, acc_local_id = accepted_name_usage_id, wfo_id = taxon_id,
    taxon_name = scientific_name, taxon_authors = scientific_name_authorship,
    taxon_status = taxonomic_status,
    ipni_id = case_when(
        wfo_id == "wfo-0001265312" ~ "545950-1",
        wfo_id == "wfo-0000605741" ~ "126223-1",
        wfo_id == "wfo-0001138844" ~ "811154-1",
        wfo_id == "wfo-0001069988" ~ "461785-1",
        wfo_id == "wfo-0001069987" ~ "461734-1",
        wfo_id == "wfo-0001098669" ~ "384680-1",
        wfo_id == "wfo-0000841098" ~ "77176713-1",
        wfo_id == "wfo-0000118782" ~ "174334-1",
        wfo_id == "wfo-0001238136" ~ "77235328-1",
        wfo_id == "wfo-0001131923" ~ "721831-1",
        wfo_id == "wfo-0000841070" ~ "288397-2",
        wfo_id == "wfo-0000989932" ~ "730628-1",
        wfo_id == "wfo-0000114078" ~ "287812-2",
        wfo_id == "wfo-0001251663" ~ "77097036-1",
        wfo_id == "wfo-0001251666" ~ "77097037-1",
        wfo_id == "wfo-0001251668" ~ "77097043-1",
        wfo_id == "wfo-0001097405" ~ "385482-1",
        wfo_id == "wfo-0001095220" ~ "686312-1",
        wfo_id == "wfo-0001214788" ~ "32126-2",
        wfo_id == "wfo-0001263203" ~ "77103256-1",
        wfo_id == "wfo-0001238147" ~ "122118-1",
        wfo_id == "wfo-0001097432" ~ "394503-1",
        wfo_id == "wfo-0001101472" ~ "693353-1",
        wfo_id == "wfo-0001253049" ~ "961931-1",
        wfo_id == "wfo-0000945834" ~ "60458350-2",
        wfo_id == "wfo-0001270351" ~ "329367-2",
        wfo_id == "wfo-0000869595" ~ "402742-1",
        wfo_id == "wfo-0000081664" ~ "208405-1",
        wfo_id == "wfo-0000700356" ~ "77309998-1",
        wfo_id == "wfo-0000703821" ~ "813065-1",
        wfo_id == "wfo-0001259250" ~ "169043-1",
        wfo_id == "wfo-0001097510" ~ "77242023-1",
        wfo_id == "wfo-0001006185" ~ "77221969-1",
        wfo_id == "wfo-0001236029" ~ "381859-1",
        wfo_id == "wfo-0001222428" ~ "647805-1",
        wfo_id == "wfo-0001088547" ~ "77243640-1",
        wfo_id == "wfo-0001115817" ~ "17259920-1",
        wfo_id == "wfo-0001010022" ~ "729042-1",
        wfo_id == "wfo-0000400261" ~ "332913-1",
        wfo_id == "wfo-0001034691" ~ "792710-1",
        wfo_id == "wfo-0001010347" ~ "735542-1",
        wfo_id == "wfo-0001006753" ~ "736860-1",
        wfo_id == "wfo-0001013208" ~ "737370-1",
        wfo_id == "wfo-0001006743" ~ "737961-1",
        wfo_id == "wfo-0001013199" ~ "738792-1",
        wfo_id == "wfo-0001115818" ~ "17030030-1",
        wfo_id == "wfo-0000935444" ~ "77239698-1",
        .default = ipni_id
      ),
    .keep = "unused"
  ) %>%
  group_by(name_to_match) %>%
  fill(ipni_id, .direction = "downup") %>%
  ungroup() %>%
  glimpse()


# UPDATE LIST OF HYBRIDS FROM WCVP BASED ON WFO DATA -----------------------------------------------

wcvp_v11_h0 %<>%
  bind_rows(
    wcvp_v11 %>%
      filter(
        name_to_match %in% wfo_h0$name_to_match | 
          (!is.na(ipni_id) & ipni_id %in% wfo_h0$ipni_id)
  ))
wcvp_v11_h0 <- wcvp_v11 %>%
  filter(plant_name_id %in% c(wcvp_v11_h0$plant_name_id, wcvp_v11_h0$accepted_plant_name_id)) %>%
  distinct()


# COMBINE HYBRID LISTS FROM WFO AND WCVP -----------------------------------------------------------

all_h0 <- wcvp_v11_h0 %>%
  transmute(
    local_id = as.character(plant_name_id),
    acc_local_id = as.character(accepted_plant_name_id), powo_id,
    ipni_id = na_if(ipni_id, "60477665-2"),
    taxon_name, taxon_authors, taxon_rank, taxon_status, name_to_match, ver, source
  ) %>%
  bind_rows(wfo_h0) %>%
  group_by(name_to_match) %>%
  fill(c(ipni_id, powo_id, wfo_id), .direction = "downup") %>%
  ungroup()


# QUERY IPNI USING NAMES ---------------------------------------------------------------------------

# add an edited full taxon name
all_h1 <- all_h0 %>%
  mutate(
    name_for_ipni = str_replace_all(taxon_name, c("× ×" = "× ", "×|\\+" = "× ")) %>%
      str_squish() %>%
      str_trim() %>%
      str_replace_all(
        c(
          "Cirsiocarduus" = "Cirsio-carduus", "Anthematricaria" = "Anthe-matricaria",
          "Anthemimatricaria" = "Anthemi-matricaria", "Epilaeliopsis" = "Epilopsis"
        )
      )
  )

# pull the names with missing IPNI IDs
names_for_ipni <- all_h1 %>%
  filter(is.na(ipni_id)) %>%
  pull(name_for_ipni)
names_for_ipni <- unique(c(names_for_ipni, str_remove_all(names_for_ipni, "-")))

# set up parallel computing
cl <- makeCluster(detectCores() - 1)
doParallel::registerDoParallel(cl)

# search IPNI
q_name_ipni_l <- foreach(
  i = seq_along(names_for_ipni),
  .packages = c("kewr"),
  .combine = "c",
  .errorhandling = "pass"
) %dopar% {
  out <- search_ipni(names_for_ipni[i])
  return(list(out))
}

# tidy the query results
q_name_ipni <- q_name_ipni_l %>%
  lapply(., kewr::tidy) %>%
  purrr::list_rbind() %>%
  clean_names() %>%
  distinct()

# match based on WFO ID
all_h1 %<>%
  left_join(
    q_name_ipni %>%
      transmute(wfo_id, ipni_id = id) %>%
      distinct(),
    by = c("wfo_id" = "wfo_id"), relationship = "many-to-many", na_matches = "never"
  ) %>%
  mutate(ipni_id = coalesce(ipni_id.x, ipni_id.y), .keep = "unused") %>%
  group_by(name_to_match) %>%
  fill(c(ipni_id, powo_id, wfo_id), .direction = "downup") %>%
  ungroup()

# match based on full names
all_h1 %<>%
  left_join(
    q_name_ipni %>%
      transmute(
        name_to_match = str_remove_all(name, "× |-") %>%
          str_c(stri_trans_general(authors, "latin-ascii"), sep = " "),
        ipni_id = id
      ) %>%
      filter(n() == 1, .by = name_to_match),
    by = c("name_to_match" = "name_to_match"),
    relationship = "many-to-many"
  ) %>%
  mutate(ipni_id = coalesce(ipni_id.x, ipni_id.y), .keep = "unused") %>%
  select(-name_for_ipni)


# QUERY IPNI USING THE IPNI IDS --------------------------------------------------------------------

# pull IPNI IDs
ipni_ids <- unique(all_h1$ipni_id)
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
write_rds(q_id_ipni_l, here("q_id_ipni_l.rds"))

# check if any IPNI IDs were skipped (shouldn't be any)
q_id_ipni_l %>%
  .[names(.) != "call"] %>%
  .[sapply(., typeof) != "list"]

# unlist IPNI query results
q_id_ipni <- q_id_ipni_l %>%
  .[sapply(., typeof) == "list"] %>%
  map(., ~ .[c("id", "wfoId", "sameCitationAs", "duplicateCitationOf")] %>%
    compact()) %>%
  data.table::rbindlist(., fill = TRUE) %>%
  mutate(ipni_id_ = id, ipni_id = id, wfo_id = wfoId, .keep = "unused")

# retrieve linked (same as, duplicate of) IPNI IDs
q_id_ipni_linked <- q_id_ipni %>%
  pivot_longer(c(sameCitationAs, duplicateCitationOf),
    values_to = "linked_ids", names_to = NULL
  ) %>%
  filter(!map_lgl(linked_ids, is.null)) %>%
  distinct() %>%
  mutate(linked_ids = map(linked_ids, ~ .[c("id", "wfoId")] %>% compact())) %>%
  unnest_wider(col = linked_ids) %>%
  rename(ipni_id_linked = id, wfo_id_linked = wfoId)

q_id_ipni_linked %<>%
  select(ipni_id_, ipni_id, ipni_id_linked) %>%
  pivot_longer(
    c(ipni_id, ipni_id_linked), names_to = NULL,
    values_drop_na = TRUE, values_to = "ipni_id"
  ) %>%
  left_join(
    q_id_ipni_linked %>%
      select(ipni_id_, wfo_id, wfo_id_linked) %>%
      pivot_longer(
        c(wfo_id, wfo_id_linked), names_to = NULL,
        values_drop_na = TRUE, values_to = "wfo_id"
      ),
    relationship = "many-to-many"
  ) %>%
  distinct()

# merge queried ID and their linked IDs
q_id_ipni %<>%
  select(ipni_id_, ipni_id, wfo_id) %>%
  filter(!(ipni_id_ %in% q_id_ipni_linked$ipni_id_)) %>%
  bind_rows(q_id_ipni_linked)

# add missing WFO taxa that are linked to the queried IPNI IDs in IPNI
setdiff(q_id_ipni$wfo_id, all_h1$local_id)
setdiff(q_id_ipni$wfo_id, wfo_v2.1$taxon_id)
setdiff(q_id_ipni$wfo_id, wfo_v3$taxon_id)
h2_add <- bind_rows(
  wfo_v2.1 %>%
    filter(taxon_id %in% setdiff(q_id_ipni$wfo_id, all_h1$local_id)),
  wfo_v3 %>%
    filter(taxon_id %in% setdiff(q_id_ipni$wfo_id, all_h1$local_id))
)
h2_add <- bind_rows(
  wfo_v2.1 %>%
    filter(taxon_id %in% c(h2_add$taxon_id, h2_add$accepted_name_usage_id)),
  wfo_v3 %>%
    filter(taxon_id %in% c(h2_add$taxon_id, h2_add$accepted_name_usage_id))
) %>%
  filter(!(taxon_id %in% all_h1$local_id)) %>%
  filter(n() == 1 | (n() > 1 & ver == "v2023-12"), .by = taxon_id)

all_h2 <- h2_add %>%
  mutate(
    local_id = taxon_id, acc_local_id = accepted_name_usage_id, wfo_id = taxon_id,
    taxon_name = scientific_name, taxon_authors = scientific_name_authorship,
    taxon_status = taxonomic_status,
    .keep = "unused"
  ) %>%
  bind_rows(all_h1)
  

# add the IPNI information to the updated WFO list
all_h2 %<>%
  left_join(
     select(q_id_ipni, ipni_id_, ipni_id) %>%
      distinct(),
    by = c("ipni_id" = "ipni_id_"), relationship = "many-to-many"
  ) %>%
  pivot_longer(c(ipni_id, ipni_id.y), values_to = "ipni_id", names_to = NULL) %>%
  distinct() %>%
  filter(n() == 1 | (n() > 1 & !is.na(ipni_id)), .by = local_id) %>%
  left_join(
      select(q_id_ipni, wfo_id, ipni_id) %>%
      distinct(),
    by = c("wfo_id" = "wfo_id"), relationship = "many-to-many", na_matches = "never"
  ) %>%
  pivot_longer(c(ipni_id.x, ipni_id.y), values_to = "ipni_id", names_to = NULL) %>%
  distinct() %>%
  filter(n() == 1 | (n() > 1 & !is.na(ipni_id)), .by = local_id) %>%
  group_by(name_to_match) %>%
  fill(c(ipni_id, powo_id, wfo_id), .direction = "downup") %>%
  ungroup() 

all_h2 <- wcvp_v11 %>%
  filter(ipni_id %in% all_h2$ipni_id & !is.na(ipni_id)) %>%
  filter(!(ipni_id %in% all_h2[all_h2$source == "WCVP", ]$ipni_id)) %>%
  transmute(
    local_id = as.character(plant_name_id),
    acc_local_id = as.character(accepted_plant_name_id), powo_id,
    ipni_id, taxon_name, taxon_authors, 
    taxon_rank, taxon_status, name_to_match, ver, source
  ) %>%
  bind_rows(all_h2) %>%
  group_by(ipni_id) %>%
  fill(c(powo_id, wfo_id), .direction = "downup") %>%
  ungroup()
  
write_rds(all_h2, here("outputs/hybrids-wfo-and-powo_0.rds"))