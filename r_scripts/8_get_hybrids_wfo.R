# Marina Golivets
# Extract hybrid names from World Flora Online (WFO)

# ==================================================================================================
# load WFO data and subset hybrids -----------------------------------------------------------------
# load the two latest versions of the WFO backbone
# I use both versions because some hybrids are excluded from the newer version
wfo_v1.1 <- data.table::fread(here("WFOTaxonomicBackbone_v.1.1_20230107/classification.txt"))
wfo_v2.1 <- data.table::fread(here("WFOTaxonomicBackbone_v.2.1_20230304/classification.csv")) %>%
  # clean raw taxon names
  mutate(scientificName = str_extract_all(scientificName, "[A-Za-z-×öëï. ]+"))
wfo_v2.1 %<>%
  unnest(scientificName) %>%
  group_by(taxonID) %>%
  summarise(scientificName = str_c(scientificName, collapse = "")) %>%
  right_join(select(wfo_v2.1, -scientificName))

setdiff(colnames(wfo_v1.1), colnames(wfo_v2.1))
setdiff(colnames(wfo_v2.1), colnames(wfo_v1.1))

# subset hybrids from the WFO backbones based on the occurrence of "×" in the names
wfo_v1.1_h0 <- wfo_v1.1 %>%
  filter(str_detect(scientificName, "×"), taxonRank != "genus")
wfo_v2.1_h0 <- wfo_v2.1 %>%
  filter(str_detect(scientificName, "×") | taxonID %in% wfo_v1.1_h0$taxonID, 
         taxonRank != "genus")
# (from v1.1 use only the names that are not in v.2.1)
wfo_v1.1_h0 %<>%
  filter(!(taxonID %in% wfo_v2.1_h0$taxonID))

wfo_h0 <- bind_rows(
  wfo_v1.1 %>%
    # subset hybrids and their accepted names (accepted names are not hybrids sometimes)
    filter(taxonID %in% c(wfo_v1.1_h0$taxonID, wfo_v1.1_h0$acceptedNameUsageID)) %>%
    mutate(wfo_version = "wfo v1.1") %>%
    janitor::clean_names(),
  wfo_v2.1 %>%
    filter(taxonID %in% c(wfo_v2.1_h0$taxonID, wfo_v2.1_h0$acceptedNameUsageID)) %>%
    mutate(wfo_version = "wfo v2.1") %>%
    janitor::clean_names()
) %>%
  group_by(taxon_id) %>%
  filter(n() == 1 | (n() > 1 & wfo_version == "wfo v2.1")) %>%
  ungroup() %>%
  mutate(
    wfo_id = taxon_id,
    ipni_id_wfo = scientific_name_id %>%
      str_remove("urn:lsid:ipni.org:names:") %>%
      if_else(str_detect(., "-"), ., NA_character_) %>%
      na_if(""),
    .keep = "unused"
  )
 
  
# fill in IPNI IDs provided in WFO (Zenodo) --------------------------------------------------------
# R.utils::gunzip(here("ipni_to_wfo.csv.gz"))
wfo_to_ipni <- data.table::fread(here("ipni_to_wfo.csv")) %>%
  mutate(ipni_id_wfo = str_remove(ipni_id, "urn:lsid:ipni.org:names:")) %>%
  group_by(wfo_id) %>%
  filter(n() == 1) %>%
  ungroup()

wfo_h0 %<>%
  left_join(wfo_to_ipni, by = c("wfo_id" = "wfo_id")) %>%
  mutate(
    ipni_id_wfo = coalesce(ipni_id_wfo.y, ipni_id_wfo.x),
    ipni_id = ipni_id_wfo,
    ipni_id_source = if_else(!is.na(ipni_id_wfo), "WFO", NA_character_),
    .keep = "unused"
  )

table(is.na(wfo_h0$ipni_id))
# FALSE  TRUE
# 19557  2098

# add IPNI IDs for some of the hybrids manually ----------------------------------------------------
wfo_h0 %<>%
  mutate(
    ipni_id_manual = case_when(
      wfo_id == "wfo-0001265312" ~ "545950-1",
      wfo_id == "wfo-0000991874" ~ "722822-1",
      wfo_id == "wfo-0000077972" ~ "77117811-1",
      wfo_id == "wfo-0000605741" ~ "126223-1",
      wfo_id == "wfo-0000423590" ~ "291314-2",
      wfo_id == "wfo-0001138844" ~ "811154-1",
      wfo_id == "wfo-0000027868" ~ "256075-1",
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
      wfo_id == "wfo-0001222424" ~ "289978-2",
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
      wfo_id == "wfo-0000816312" ~"290198-2",
      wfo_id == "wfo-0000954861" ~ "592901-1",
      wfo_id == "wfo-0000178123" ~ "144733-2",
      wfo_id == "wfo-0001056138" ~ "144733-2",
      wfo_id == "wfo-0001115818" ~ "17030030-1",
      wfo_id == "wfo-0000096024" ~ "287578-2",
      "wfo-0000693108" ~ "672862-1",
      "wfo-0000935444" ~ "77239698-1",
      "wfo-0001107413" ~ "324133-2",
      
    ),
    ipni_id = coalesce(ipni_id_manual, ipni_id),
    ipni_id_source = if_else(
      !is.na(ipni_id_manual), "manual search", ipni_id_source
    )
  )

# remove duplicates
wfo_h0 %<>%
  group_by(ipni_id) %>%
  filter(is.na(ipni_id) | n() == 1 | n() > 1 & wfo_version == "wfo v2.1") %>%
  filter(is.na(ipni_id) | n() == 1 | n() > 1 & taxonomic_status == "Accepted") %>%
  ungroup()
table(wfo_h0$ipni_id_source)
# manual search   WFO 
# 28              19488 

table(is.na(wfo_h0$ipni_id))
# FALSE  TRUE 
# 19516  2053


# ==================================================================================================
# query IPNI using the IDs available in WFO --------------------------------------------------------
# set up parallel computing
cl <- makeCluster(detectCores() - 1)
doParallel::registerDoParallel(cl)

# query IPNI using IPNI ID as provided in WFO
ipni_id <- unique(wfo_h0$ipni_id[!is.na(wfo_h0$ipni_id)])
q_id_ipni_l <- foreach(
  i = seq_along(ipni_id),
  .packages = c("kewr"),
  .combine = "c",
  .errorhandling = "pass"
) %dopar% {
  out <- lookup_ipni(ipni_id[i], type = "taxon")
  return(list(out))
}

# check which IPNI IDs were skipped
q_id_ipni_l %>%
  .[names(.) != "call"] %>%
  .[sapply(., typeof) != "list"]

# unlist IPNI query results
q_id_ipni <- q_id_ipni_l %>%
  .[sapply(., typeof) == "list"] %>%
  map(., ~ .[c("suppressed", "inPowo", "topCopy", "id", "wfoId")]) %>%
  data.table::rbindlist(., fill = TRUE) %>%
  clean_names() %>%
  select(-na)

# add missing WFO taxa that are linked to the queried IPNI IDs in IPNI
wfo_h0_add <- wfo_v2.1 %>%
  filter(taxonID %in% setdiff(q_id_ipni$wfo_id, wfo_h0$wfo_id)) 
wfo_h0_add <- wfo_v2.1 %>%
  filter(taxonID %in% c(wfo_h0_add$taxonID, wfo_h0_add$acceptedNameUsageID)) %>%
  filter(!(taxonID %in% wfo_h0$wfo_id))

wfo_h0_ <- wfo_h0_add %>%
  janitor::clean_names() %>%
  mutate(
    wfo_id = taxon_id,
    ipni_id_wfo = scientific_name_id %>%
      str_remove("urn:lsid:ipni.org:names:") %>%
      na_if(""),
    wfo_version = "wfo v2.1",
    .keep = "unused"
  ) %>%
  bind_rows(wfo_h0) %>%
  group_by(wfo_id) %>%
  filter(n() == 1 | (n() > 1 & wfo_version == "wfo v2.1")) %>%
  ungroup()

# add the IPNI information to the updated WFO list
wfo_h1 <- q_id_ipni %>%
  transmute(wfo_id, ipni_id_ipni = id) %>%
  right_join(wfo_h0_) %>%
  mutate(
    ipni_id_source = if_else(is.na(ipni_id) & !is.na(ipni_id_ipni), "IPNI", ipni_id_source),
    ipni_id = coalesce(ipni_id, ipni_id_ipni)
  ) %>%
  left_join(
    q_id_ipni %>%
      select(-wfo_id) %>%
      rename(ipni_id = id)
  )

# remove suppressed IDs and IDs that are not top copy (to be replaced with new IDs)
wfo_h1 %<>%
  mutate(
    ipni_id = case_when(
      suppressed == TRUE & in_powo == FALSE ~ NA_character_,
      top_copy == FALSE & in_powo == FALSE ~ NA_character_,
      .default = ipni_id
    ),
    across(c(top_copy, suppressed, in_powo, ipni_id_source), ~ if_else(is.na(ipni_id), NA, .))
  )

table(is.na(wfo_h1$ipni_id))
# FALSE  TRUE 
# 19456  2344

# query IPNI using WFO names -----------------------------------------------------------------------
wfo_h2 <- wfo_h1 %>%
  transmute(
    wfo_id, ipni_id, top_copy, suppressed, in_powo, ipni_id_source,
    wfo_name_v = scientific_name,
    wfo_name = scientific_name %>%
      str_replace_all(c("× ×" = "× ", "×" = "× ")) %>%
      str_squish() %>%
      str_trim(),
    wfo_authors_v = str_replace_all(
      scientific_name_authorship,
      c(
        "P‚nzes" = "Pénzes",
        "Palm,n" = "Palmén",
        "Andr‚" = "André",
        "„" = "ä",
        "H,rincq" = "Hérincq",
        "Kr?ssm" = "Krüssm",
        "¢" = "ó",
        '"' = "ö"
      )
    ),
    wfo_rank = str_to_lower(taxon_rank),
    wfo_genus = genus,
    wfo_family = family,
    wfo_status = str_to_lower(taxonomic_status),
    wfo_acc_id = if_else(wfo_status == "accepted", wfo_id, accepted_name_usage_id),
    wfo_version
  ) %>%
  # manually edit some of the genera names
  mutate(
    across(c(wfo_name, wfo_genus), ~ str_replace_all(
      .,
      c(
        # "-" = "",
        "Cirsiocarduus" = "Cirsio-carduus",
        "Anthematricaria" = "Anthe-matricaria",
        "Anthemimatricaria" = "Anthemi-matricaria",
        "Crepi-Hieracium" = "Crepi-hieracium",
        "Loroglorchis" = "Lorogl-orchis",
        "Crataegomespilus" = "Crataego-mespilus",
        "Epilaeliopsis" = "Epilopsis"
      )
    )),
    wfo_authors = stringi::stri_trans_general(wfo_authors_v, "latin-ascii") %>%
      str_replace_all(
        c(
          "Hal csy" = "Halacsy",
          "G yer" = "Gayer",
          "Borb s" = "Borbas",
          "Prod n" = "Prodan",
          "Soj k" = "Sojak",
          "J v" = "Jav",
          ",p\\.p.*" = "",
          ",$" = ""
        )
      ),
    # replace blanks with NAs
    across(where(is.character), ~ na_if(., ""))
  )

# query IPNI using canonical names (only for names that are missing IPNI IDs)
ipni_name <- unique(
  c(wfo_h2$wfo_name[is.na(wfo_h2$ipni_id)],
    wfo_h2$wfo_name[is.na(wfo_h2$ipni_id)] %>% str_remove_all("-")
    ))
q_name_ipni_l <- foreach(
  i = seq_along(ipni_name),
  .packages = c("kewr"),
  .combine = "c",
  .errorhandling = "pass"
) %dopar% {
  out <- search_ipni(ipni_name[i])
  return(list(out))
}

q_name_ipni <- q_name_ipni_l %>%
  lapply(., kewr::tidy) %>%
  purrr::list_rbind() %>%
  janitor::clean_names()

# match based on WFO ID
wfo_h2 %<>%
  left_join(
    q_name_ipni %>%
      transmute(
        wfo_id, top_copy, in_powo, suppressed, id,
        ipni_id_source = "IPNI"
        ),
    by = c("wfo_id" = "wfo_id")
  ) %>%
  mutate(
    ipni_id = coalesce(ipni_id, id),
    top_copy = coalesce(top_copy.x, top_copy.y),
    in_powo = coalesce(in_powo.x, in_powo.y),
    suppressed = coalesce(suppressed.x, suppressed.y),
    ipni_id_source = coalesce(ipni_id_source.x, ipni_id_source.y),
    .keep = "unused"
  ) %>%
  distinct() %>%
  group_by(wfo_id) %>%
  filter(n() == 1 | (n() > 1 & in_powo == TRUE)) %>%
  ungroup()

table(is.na(wfo_h2$ipni_id))
# FALSE  TRUE
# 19530  2270  

# match based on full names
wfo_h2 %<>%
  mutate(
    name_to_match = str_remove_all(wfo_name, "×|-") %>%
    str_squish() %>%
    str_trim()
    ) %>%
  left_join(
    q_name_ipni %>%
      transmute(
        id, top_copy, in_powo, suppressed,
        name_to_match = str_remove_all(name, "× |-"),
        authors_to_match = stringi::stri_trans_general(authors, "latin-ascii"),
        ipni_id_source = "IPNI"
      ) %>%
      group_by(name_to_match, authors_to_match) %>%
      filter(n() == 1),
    by = c("name_to_match" = "name_to_match", "wfo_authors" = "authors_to_match")
  ) %>%
  mutate(
    ipni_id = coalesce(ipni_id, id),
    top_copy = coalesce(top_copy.x, top_copy.y),
    in_powo = coalesce(in_powo.x, in_powo.y),
    suppressed = coalesce(suppressed.x, suppressed.y),
    ipni_id_source = coalesce(ipni_id_source.x, ipni_id_source.y),
    .keep = "unused"
  )

table(is.na(wfo_h2$ipni_id))
# FALSE  TRUE
# 20242  1558  

# match based on canonical names
wfo_h2 %<>%
  left_join(
    q_name_ipni %>%
      transmute(
        id, top_copy, in_powo, suppressed,
        name_to_match = str_remove_all(name, "× |-"),
        ipni_id_source = "IPNI"
      ) %>%
      group_by(name_to_match) %>%
      filter(n() == 1),
    by = c("name_to_match" = "name_to_match")
  ) %>%
  mutate(
    ipni_id = coalesce(ipni_id, id),
    top_copy = coalesce(top_copy.x, top_copy.y),
    in_powo = coalesce(in_powo.x, in_powo.y),
    suppressed = coalesce(suppressed.x, suppressed.y),
    ipni_id_source = coalesce(ipni_id_source.x, ipni_id_source.y),
    .keep = "unused"
  )
table(is.na(wfo_h2$ipni_id))
# FALSE  TRUE
# 20597  1203


# ==================================================================================================
# query POWO for taxa with unmatched IPNI IDs using canonical names --------------------------------
wfo_wo_ipni <- unique(
  c(
    wfo_h2$wfo_name[is.na(wfo_h2$ipni_id)],
    wfo_h2$wfo_name[is.na(wfo_h2$ipni_id)] %>% str_remove_all("-")
  ))

q_name_powo_l <- foreach(
  i = seq_along(wfo_wo_ipni),
  .packages = c("taxize"),
  .combine = "c",
  .errorhandling = "pass"
) %dopar% {
  Sys.sleep(.2)
  out <- get_pow_(wfo_wo_ipni[i], messages = FALSE)
  return(out)
}

q_name_powo <- q_name_powo_l %>%
  .[sapply(., typeof) == "list"] %>%
  bind_rows(.id = "wfo_name") %>%
  select(-images) %>%
  distinct() %>%
  janitor::clean_names()

# match WFO to POWO using full names
wfo_h3 <- wfo_h2 %>%
  mutate(name_pow = paste(str_remove_all(wfo_name, "× |-"), wfo_authors)) %>%
  left_join(
    q_name_powo %>%
      transmute(
        ipni_id_pow = str_remove(fq_id, "urn:lsid:ipni.org:names:"),
        name_pow = paste(
          str_remove_all(name, "× |-"),
          stringi::stri_trans_general(author, "latin-ascii")
        ),
        ipni_id_source = "POWO"
      ) %>%
      group_by(name_pow) %>%
      filter(n() == 1),
    by = c("name_pow" = "name_pow")
  ) %>%
  mutate(ipni_id = coalesce(ipni_id, ipni_id_pow)) %>% 
  mutate(ipni_id_source = coalesce(ipni_id_source.x, ipni_id_source.y),
         .keep = "unused") %>%
  select(-name_pow) %>%
  ungroup()
table(is.na(wfo_h3$ipni_id))
# FALSE  TRUE
# 21363   437 

# match WFO to POWO using canonical names
wfo_h3 %<>%
  mutate(name_pow = str_remove_all(wfo_name, "× ")) %>%
  left_join(
    q_name_powo %>%
      transmute(
        ipni_id_pow = str_remove(fq_id, "urn:lsid:ipni.org:names:"),
        name_pow = str_remove_all(name, "× |-"),
        ipni_id_source = "POWO"
      ) %>%
      group_by(name_pow) %>%
      filter(n() == 1),
    by = c("name_pow" = "name_pow")
  ) %>%
  mutate(
    ipni_id_pow = coalesce(ipni_id_pow.x, ipni_id_pow.y),
    ipni_id = coalesce(ipni_id, ipni_id_pow),
         ipni_id_source = coalesce(ipni_id_source.x, ipni_id_source.y),
         .keep = "unused"
  )
table(is.na(wfo_h3$ipni_id))
# FALSE  TRUE
# 21478   322 



