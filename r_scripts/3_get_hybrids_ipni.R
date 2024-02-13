# Marina Golivets
# BiCIKL Hackathon (Sept 20-24, 2021, Meise, Belgium)
# Extract hybrid and parent taxon names from IPNI


# get hybrid names from IPNI -----------------------------------------------------------------------
hybr_ipni0 <- data.table::fread(here("input/hybrids_ipni.csv"), na.strings = c(""))
glimpse(hybr_ipni0)
unique(hybr_ipni0$rank_s_alphanum)

hybr_ipni <- hybr_ipni0 %>%
  transmute(
    ipni_id = id,
    rank_s_alphanum, hybrid_parents_s_lower, lookup_hybrid_parent_id,
    verbatim = paste(taxon_scientific_name_s_lower, authors_t) %>%
      str_squish() %>%
      str_trim(),
    type = "name"
  ) %>%
  # remove taxa higher than species
  filter(!(rank_s_alphanum %in% c(
    "gen.", "nothosect.", "nothoser.", "nothosubgen.",
    "nothosubsect.", "sect.", "subsect."
  )))
hybr_ipni_excl <- hybr_ipni0 %>%
  filter(!(id %in% hybr_ipni$ipni_id)) %>%
  mutate(reason = "higher rank")

hybr_ipni <- bind_cols(
  hybr_ipni %>%
    mutate(
      parent = str_remove_all(hybrid_parents_s_lower, pattern = "^\\(|\\)$") %>%
        str_squish() %>%
        str_trim() %>%
        str_split(" × (?=[A-Z])"), .keep = "unused"
    ) %>%
    unnest(cols = "parent"),
  hybr_ipni %>%
    mutate(parent_ipni_id = str_split(lookup_hybrid_parent_id, ","), .keep = "unused") %>%
    unnest("parent_ipni_id") %>%
    select(parent_ipni_id)
) %>%
  distinct() %>%
  mutate(
    parent_ipni_id = if_else(
      is.na(parent), NA_character_, paste0("urn:lsid:ipni.org:names:", parent_ipni_id)
    ),
    parent = str_remove(parent, "^× ") %>%
      case_match(
        "Helianthemum almeriense Pau subsp. almeriense" ~ "Helianthemum almeriense Pau",
        "Ophrys sipotensis" ~ "Ophrys sipontensis",
        "Sideritis ibanyezii" ~ "Sideritis ibanezii",
        "Cistus lauriflorus" ~ "Cistus laurifolius",
        "Ophrys oxyrrhynchos subsp. oxyrrhynchos" ~ "Ophrys oxyrrhynchos",
        "Globulea nudicaulis Haw. var. nudicaulis" ~ "Globulea nudicaulis",
        "Eremaea asterocarpa Hnatiuk var. asterocarpa" ~ "Eremaea asterocarpa Hnatiuk",
        "Ophrys bertoloniiformis O.Danesch & E.Danesch subsp. bertoloniiformis" ~ "Ophrys bertoloniiformis O.Danesch & E.Danesch",
        "Potentilla josphiana" ~ "Potentilla josephiana",
        "Septas capensis L. var. capensis" ~ "Septas capensis L.",
        "Crataegus ambigiua subsp. ambigua" ~ "Crataegus ambigiua",
        "Rhododendron kaepferi var. macrogemma" ~ "Rhododendron kaempferi var. macrogemmum",
        "Cattleya leopoldi" ~ "Cattleya leopoldii",
        "Galactites durieui" ~ "Galactites duriaei Spach",
        "Monanthes laxiflora Bolle var. laxiflora" ~ "Monanthes laxiflora Bolle",
        "Baccharis santelicis Phil. subsp. santelicis" ~ "Baccharis santelicis Phil.",
        "Baccharis elatoides" ~ "Baccharis elaeoides",
        "Saxifraga continentales" ~ "Saxifraga continentalis",
        "Glyceria ischyromeura" ~ "Glyceria ischyroneura",
        "Ramonda myconii" ~ "Ramonda myconi",
        "Amygdalus elaeagnifolia subsp. elaeagnifolia" ~ "Prunus elaeagrifolia",
        "Salix salvifolia" ~ "Salix salviifolia",
        "Cardamine matthiolii" ~ "Cardamine matthioli",
        "Aeonium spathulatum (Hornem.) Praeger var. spathulatum" ~ "Aeonium spathulatum (Hornem.) Praeger",
        "Aeonium palmensis" ~ "Aeonium palmense",
        "Achillea bibersteinii" ~ "Achillea biebersteinii",
        "Cistus symphytifolius Lam. var. symphytifolius" ~ "Cistus symphytifolius Lam.",
        "Reynoutria sachaliensis" ~ "Reynoutria sachalinensis",
        "Ophrys oestrifera M.Bieb. subsp. oestrifera" ~ "Ophrys oestrifera M.Bieb.",
        "Gentiana campestris Geners. subsp. campestris" ~ "Gentiana campestris Geners.",
        "Cerasus fructiosa" ~ "Cerasus fruticosa",
        "Rhinanthus alpinus Lam. subsp. alpinus" ~ "Rhinanthus alpinus Lam.",
        "Avenula gervaisii Holub subsp. gervaisii" ~ "Avenula gervaisii Holub",
        "Orchis anatolica Boiss. subsp. anatolica" ~ "Orchis anatolica Boiss.",
        "Nepenthes thorelli" ~ "Nepenthes thorelii",
        "Solanum acule" ~ "Solanum aculeastrum",
        "Ophrys holoserica (Burm.f.) Greuter subsp. holosericea" ~ "Ophrys holoserica (Burm.f.) Greuter",
        "Orchis obliensis" ~ "Orchis olbiensis",
        "Halophila nipponica J.Kuo subsp. nipponica" ~ "Halophila nipponica J.Kuo",
        "Daphne cneorum L. subsp. cneorum" ~ "Daphne cneorum L.",
        "Aeonium canariense Webb & Berthel. var. canariense" ~ "Aeonium canariense Webb & Berthel.",
        "Aeonium urbicum (C.Sm. ex Hornem.) Webb & Berthel. var. urbicum" ~ "Aeonium urbicum (C.Sm. ex Hornem.) Webb & Berthel.",
        "Pinus densata Masters subsp. densata" ~ "Pinus densata Masters",
        "Salicornia persica Akhani subsp. persica" ~ "Salicornia persica Akhani",
        "Satureja obovata Briq. subsp. obovata" ~ "Satureja obovata Briq.",
        "Potentilla stenomalla" ~ "Potentilla stenophylla",
        "Ophrys ciliata Biv. subsp. ciliata" ~ "Ophrys ciliata Biv.",
        "Agrostis castellana Boiss. & Reut. var. castellana" ~ "Agrostis castellana Boiss.",
        "Woodsia oregana D.C.Eaton f. oregana T.M.C.Taylor" ~ "Woodsia oregana D.C.Eaton",
        "Cerasus apetala (Siebold & Zucc.) Ohle var. apetala" ~ "Cerasus apetala (Siebold & Zucc.) Ohle",
        "Cerasus incisa Loisel. var. incisa" ~ "Cerasus incisa Loisel.",
        "Rhododendron smirnovii" ~ "Rhododendron smirnowii",
        "Viola nemoralis Kuetz. subsp. nemoralis" ~ "Viola nemoralis Kuetz.",
        "Pachypodium densiflorum Baker subsp. densiflorum" ~ "Pachypodium densiflorum Baker",
        "Rubus carpinifoluis" ~ "Rubus carpinifolius",
        "Pseudolysimachion incanum (L.) Holub subsp. incanum" ~ "Pseudolysimachion incanum (L.) Holub",
        "Citrus limetta Risso subsp. limetta" ~ "Citrus limetta Risso",
        "Globularia repens subsp. repens" ~ "Globularia repens",
        "Crataegus ambigiua" ~ "Crataegus ambigua",
        .default = parent
      ),
    edited = verbatim %>%
      str_remove("pro sp.") %>%
      str_squish() %>%
      str_trim(),
    edited2 = str_replace_all(verbatim, c(" × " = " ×", "^× " = "×")),
    source = "ipni"
  )

# retrieve data from POW ---------------------------------------------------------------------------
pow_h_l <- pbapply::pbsapply(
  unique(hybr_ipni$ipni_id),
  function(x) {
    Sys.sleep(.5)
    try(taxize::pow_lookup(x, include = c("distribution", "descriptions")))
  },
  cl = cl, USE.NAMES = FALSE
)

pow_h <- pow_h_l[sapply(pow_h_l, length) == 2] %>%
  lapply(., function(x) x$meta) %>%
  lapply(., rlist::list.flatten) %>%
  lapply(., function(x) x[sapply(x, length) == 1]) %>%
  data.table::rbindlist(., fill = TRUE) %>%
  janitor::clean_names() %>%
  distinct()

pow_h_l_add <- pbapply::pbsapply(
  setdiff(unique(hybr_ipni$ipni_id), pow_h$fq_id),
  function(x) {
    Sys.sleep(.5)
    try(taxize::pow_lookup(x, include = c("distribution", "descriptions")))
  },
  cl = cl, USE.NAMES = FALSE
)

pow_h_add <- pow_h_l_add[sapply(pow_h_l_add, length) == 2] %>%
  lapply(., function(x) x$meta) %>%
  lapply(., rlist::list.flatten) %>%
  lapply(., function(x) x[sapply(x, length) == 1]) %>%
  data.table::rbindlist(., fill = TRUE) %>%
  janitor::clean_names() %>%
  distinct()

pow_h %<>% bind_rows(pow_h_add)

# exclude artificial hybrids
hybr_ipni_excl %<>%
  bind_rows(
  hybr_ipni0 %>%
  filter(id %in% (pow_h %>% 
             filter(taxonomic_status == "Artificial_Hybrid") %>%
             pull(fq_id))) %>%
  mutate(reason = "artificial hybrid")
  )
hybr_ipni %<>% filter(!(ipni_id %in% hybr_ipni_excl$id))

tmp <- hybr_ipni %>%
  filter(!(ipni_id %in% (pow_h %>% pull(fq_id))))
table(tmp$rank_s_alphanum)

"Streptocarpus denticulatus S. dunnii" 
"Symphytum bulbosum S. tuberosum"
"Hemionitis acrostica H. pteridioides"
"Huperzia lucidula H. porophila"  
"Carduus acanthoides C. nigrescens"  
"Salix alba S. nigra" 
"Platanthera leucophaea P. psycodes" 
"Sticherus ferrugineus S. furcatus"
"Asplenium protovirchowii A. virchowii"    

tmp <- pow_h %>%
  mutate(
  parent = str_remove_all(hybrid_formula, pattern = "^\\(|\\)$") %>%
    str_split(pattern = "×| x |\\+")
) %>%
  unnest(cols = "parent") %>%
mutate(
  parent = parent %>%
    str_squish() %>% 
    str_trim(),
  parent = case_when(
    grepl("perhaps|\\?| sp\\.$|derived", parent) ~ NA_character_,
    grepl("^[a-z]", parent) ~ paste(genus, parent),
    grepl("^[A-Z]\\.|^[A-Z] ", parent) ~ paste(genus, gsub("^[A-Z]\\.||^[A-Z]", "", parent)),
    TRUE ~ parent
  ) %>%
    str_squish() %>% 
    str_trim())

gbif_p <- gbif.cmpfn(
  unique(tmp$parent),
  seq_along(unique(tmp$parent))
) %>%
  filter(str_detect(scientificName, "×| x ") | confidence > 95 | confidence == 0) %>%
  group_by(name) %>%
  filter(n() == 1) %>%
  ungroup()
setdiff(unique(tmp$parent), gbif_p$name)

# filter(
#   str_detect(edited, "'| hort\\.|Hort\\.", negate = TRUE),  # exclude horticultural varieties
# !(parent %in% c(
#   "Not given", "Gymnocalycium ?", "Potentilla stenomalla", "Phyllocactus sp",
#   "Saxifraga s", "Tilia t", "Gr. nivalis", "Or uniflora", "Cf. panizzianus",
#   "O. aphroditae"
# ))
# )

# match names to GBIF ------------------------------------------------------------------------------
gbif_h <- gbif.cmpfn(
  unique(hybr_ipni$edited2),
  seq_along(unique(hybr_ipni$edited2))
) %>%
  filter(str_detect(scientificName, "×| x ") | confidence > 95 | confidence == 0) %>%
  group_by(name) %>%
  filter(n() == 1) %>%
  ungroup()
setdiff(unique(hybr_ipni$edited2), gbif_h$name)

hybr_ipni %<>%
  left_join(
    gbif_h %>%
      transmute(
        name, usageKey,
        acceptedUsageKey = coalesce(acceptedUsageKey, usageKey)
      ),
    by = c("edited2" = "name")
  )

gbif_p <- gbif.cmpfn(
  unique(hybr_ipni$parent),
  seq_along(unique(hybr_ipni$parent))
) %>%
  filter(str_detect(scientificName, "×| x ") | confidence > 95 | confidence == 0) %>%
  group_by(name) %>%
  filter(n() == 1) %>%
  ungroup()
setdiff(unique(hybr_ipni$parent), gbif_p$name)

hybr_ipni_2p <- hybr_ipni %>%
  group_by(edited) %>%
  filter(n() == 2) %>%
  select(-parent_taxonID) %>%
  mutate(parent_type = c("parent1", "parent2")) %>%
  ungroup() %>%
  pivot_wider(values_from = parent, names_from = parent_type)

hybr_ipni_2p %<>%
  left_join(
    gbif_p %>%
      transmute(
        parent1 = name,
        parent1_usageKey = usageKey,
        parent1_acceptedUsageKey = coalesce(acceptedUsageKey, usageKey)
      )
  ) %>%
  left_join(
    gbif_p %>%
      transmute(
        parent2 = name,
        parent2_usageKey = usageKey,
        parent2_acceptedUsageKey = coalesce(acceptedUsageKey, usageKey)
      )
  ) %>%
  distinct()

hybr_ipni_2p %<>%
  bind_rows(
    hybr_ipni %>%
      group_by(edited) %>%
      filter(n() == 1) %>%
      filter(is.na(parent))
  )




write_csv2(here("output/hybrids_ipni.csv"))
