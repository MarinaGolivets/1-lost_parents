


# query IPNI using names from WFO ------------------------------------------------------------------
# use family names in WFO to query IPNI
f_ipni_l <- unique(wfo_h$family) %>%
  pbapply::pbsapply(
    .,
    function(x) {
      Sys.sleep(.5)
      try(taxize::ipni_search(family = x, output = "extended"))
    },
    cl = cl, USE.NAMES = FALSE
  )

# use genus names in WFO to query IPNI
g_ipni_l <- unique(wfo_h$genus_edited) %>%
  pbapply::pbsapply(
    .,
    function(x) {
      Sys.sleep(.5)
      try(taxize::ipni_search(genus = x, output = "extended"))
    },
    cl = cl, USE.NAMES = FALSE
  )

# combine the query results and store them as a tibble
ipni <- c(f_ipni_l, g_ipni_l) %>%
  lapply(., function(x) x[sapply(x, length) != 0]) %>%
  data.table::rbindlist(., fill = TRUE) %>%
  select(
    id, full_name_without_family_and_authors, authors,
    original_hybrid_parentage, original_remarks, rank
  ) %>%
  mutate(across(everything(), ~ na_if(., ""))) %>%
  distinct()

# match WFO to IPNI using full names
wfo_h %<>%
  mutate(full_name_ipni = paste(str_replace_all(name_edited, "×", "× "), authors)) %>%
  left_join(
    ipni %>%
      transmute(
        id_ = id,
        full_name_ipni = paste(full_name_without_family_and_authors, authors)
      ) %>%
      group_by(full_name_ipni) %>%
      filter(n() == 1)
  ) %>%
  mutate(ipni_id = coalesce(ipni_id, id_), .keep = "unused") %>%
  distinct()
table(is.na(wfo_h$ipni_id))
# FALSE  TRUE
# 15746  4701

# wfo_h %<>%
#   mutate(full_name_ipni = paste(str_remove_all(name_edited, "×"), authors)) %>%
#   left_join(
#     ipni %>%
#       transmute(
#       id_ = id,
#       full_name_ipni = paste(
#         str_remove_all(full_name_without_family_and_authors, "× ")
#       ),
#       rank
#     ) %>%
#       group_by(full_name_ipni) %>%
#       filter(n() == 1)
#   ) %>%
#   mutate(ipni_id = coalesce(ipni_id, id_), .keep = "unused")
# table(is.na(wfo_h$ipni_id))
# FALSE  TRUE
# 15743  4704

# match WFO to IPNI using canonical names
wfo_h %<>%
  mutate(name_ipni = str_remove_all(name_edited, "×")) %>%
  left_join(
    ipni %>%
      transmute(
        id_ = id,
        name_ipni = str_remove_all(full_name_without_family_and_authors, "× ")
      ) %>%
      group_by(name_ipni) %>%
      filter(n() == 1)
  ) %>%
  mutate(ipni_id = coalesce(ipni_id, id_), .keep = "unused")
table(is.na(wfo_h$ipni_id))
# FALSE  TRUE
# 18276  2171


# add hybrid names from IPNI that are not in WFO
all_h <- bind_rows(
  wfo_h,
  ipni %>%
    filter(
      str_detect(full_name_without_family_and_authors, "× "),
      !(id %in% wfo_h$ipni_id),
      str_count(full_name_without_family_and_authors, "\\S+") > 2
    ) %>%
    transmute(
      id,
      ipni_id = id,
      name = full_name_without_family_and_authors,
      name_edited = str_replace_all(full_name_without_family_and_authors, "× ", "×"),
      authors, rank,
      source = "ipni"
    ) %>%
    distinct()
) %>%
  left_join(
    ipni %>%
      transmute(
        ipni_id = id,
        formula_ipni = original_hybrid_parentage,
        remarks_ipni = original_remarks
      ) %>%
      distinct()
  ) %>%
  mutate(
    ipni_id = if_else(
      !is.na(ipni_id),
      paste0("urn:lsid:ipni.org:names:", ipni_id), NA_character_
    )
  )

# note that there are some duplicated names in the list of hybrids because IPNI provides multiple
# IDs for the same species sometimes

# filter(!(ipni_id %in% c(
#   "2928272-1", "3013719-1", "2928273-1", "249825-1", "639579-1", "999831-1",
#   "359888-1", "3000278-1", "68622-1", "276437-1", "60477048-2", "60477047-2",
#   "56798-1", "60461580-2", "77192251-1", "77071625-1", "940138-1", "77222609-1",
#   "2625853-1", "2990725-1", "2962183-4", "77159367-1", "287676-2", "741957-1"
# )))
# 77248466-1
# 77248334-1


# add IPNI IDs from the data provided at the hackathon ---------------------------------------------
# load the IPNI hybrids subset
ipni_h <- data.table::fread(here("input/hybrids_ipni.csv"), na.strings = c("")) %>%
  filter(str_detect(id, "212166-1", negate = TRUE)) %>%
  # remove taxa higher than species
  filter(!(rank_s_alphanum %in% c(
    "gen.", "nothosect.", "nothoser.", "nothosubgen.",
    "nothosubsect.", "sect.", "subsect."
  )))

# match WFO to IPNI using full names
all_h %<>%
  mutate(full_name_ipni = paste(str_replace_all(name, "×", "× "), authors)) %>%
  left_join(
    ipni_h %>%
      transmute(
        id_ = id,
        full_name_ipni = paste(taxon_scientific_name_s_lower, authors_t)
      ) %>%
      group_by(full_name_ipni) %>%
      filter(n() == 1)
  ) %>%
  mutate(ipni_id = coalesce(ipni_id, id_), .keep = "unused") %>%
  distinct()
table(is.na(all_h$ipni_id))
# FALSE  TRUE
# 23711  1999

all_h %<>%
  mutate(full_name_ipni = paste(str_remove_all(name_edited, "×"), authors)) %>%
  left_join(
    ipni_h %>%
      transmute(
        id_ = id,
        full_name_ipni = paste(
          str_remove_all(taxon_scientific_name_s_lower, "× "),
          authors_t
        )
      ) %>%
      group_by(full_name_ipni) %>%
      filter(n() == 1)
  ) %>%
  mutate(ipni_id = coalesce(ipni_id, id_), .keep = "unused")
table(is.na(all_h$ipni_id))
# FALSE  TRUE
# 23712  1998

# match WFO to IPNI using canonical names
all_h %<>%
  mutate(name_ipni = str_remove_all(name_edited, "×")) %>%
  left_join(
    ipni_h %>%
      transmute(
        id_ = id,
        name_ipni = str_remove_all(taxon_scientific_name_s_lower, "× ")
      ) %>%
      group_by(name_ipni) %>%
      filter(n() == 1)
  ) %>%
  mutate(ipni_id = coalesce(ipni_id, id_), .keep = "unused")
table(is.na(all_h$ipni_id))
# FALSE  TRUE
# 23857  1853

# add hybrid names from IPNI that are not in WFO
all_h %<>%
  bind_rows(
    ipni_h %>%
      filter(
        str_detect(taxon_scientific_name_s_lower, "× "),
        !(id %in% all_h$ipni_id)
      ) %>%
      transmute(
        id,
        ipni_id = id,
        name = taxon_scientific_name_s_lower,
        name_edited = str_replace_all(taxon_scientific_name_s_lower, "× ", "×"),
        authors = authors_t,
        rank = rank_s_alphanum,
        source = "ipni"
      ) %>%
      distinct()
  ) %>%
  left_join(
    ipni_h %>%
      transmute(
        ipni_id = id,
        parents_ipni = hybrid_parents_s_lower
      ) %>%
      distinct()
  )
table(is.na(all_h$ipni_id))
# FALSE  TRUE
# 24327  1853

# get POW taxon names and IDs based on canonical name matches --------------------------------------
pow_l <- all_h %>%
  filter(is.na(ipni_id)) %>%
  pull(name) %>%
  str_replace_all("×", "× ") %>%
  str_squish() %>%
  pbapply::pbsapply(
    .,
    function(x) {
      Sys.sleep(.5)
      try(taxize::get_pow_(x, messages = FALSE))
    },
    cl = cl, USE.NAMES = FALSE
  )

pow_l_add <- all_h %>%
  filter(is.na(ipni_id)) %>%
  pull(name) %>%
  str_replace_all("×", "× ") %>%
  str_squish() %>%
  .[sapply(pow_l, length) == 1] %>%
  pbapply::pbsapply(
    .,
    function(x) {
      Sys.sleep(.5)
      try(taxize::get_pow_(x, messages = FALSE))
    },
    cl = cl, USE.NAMES = FALSE
  )

pow <- c(pow_l, pow_l_add) %>%
  lapply(., function(x) x[sapply(x, length) > 0]) %>%
  .[sapply(., length) > 1] %>%
  bind_rows() %>%
  select(-images) %>%
  distinct()

# match WFO to POWO using full names
all_h %<>%
  mutate(name_pow = paste(str_remove_all(name, "×"), authors)) %>%
  left_join(
    pow %>%
      transmute(
        fqId,
        name_pow = paste(str_remove_all(name, "× "), author)
      ) %>%
      group_by(name_pow) %>%
      filter(n() == 1)
  ) %>%
  mutate(ipni_id = coalesce(ipni_id, fqId), .keep = "unused")
table(is.na(all_h$ipni_id))
# FALSE  TRUE
# 25627   553

# match WFO to POWO using canonical names
all_h %<>%
  mutate(name_pow = str_remove_all(name, "×")) %>%
  left_join(
    pow %>%
      transmute(
        fqId,
        name_pow = str_remove_all(name, "× ")
      ) %>%
      group_by(name_pow) %>%
      filter(n() == 1)
  ) %>%
  mutate(ipni_id = coalesce(ipni_id, fqId), .keep = "unused")
table(is.na(all_h$ipni_id))
# FALSE  TRUE
# 25835   345

# retrieve information from POWO -------------------------------------------------------------------
pow_info_l <- pbapply::pbsapply(
  unique(all_h$ipni_id),
  function(x) {
    Sys.sleep(.2)
    try(taxize::pow_lookup(x, include = c("distribution")))
  },
  cl = cl, USE.NAMES = FALSE
)

all_h_pow <- pow_info_l[sapply(pow_info_l, length) == 2] %>%
  lapply(., function(x) x$meta) %>%
  lapply(., rlist::list.flatten) %>%
  lapply(., function(x) x[sapply(x, length) == 1]) %>%
  data.table::rbindlist(., fill = TRUE) %>%
  janitor::clean_names() %>%
  distinct()

pow_info_l_add <- pbapply::pbsapply(
  setdiff(unique(all_h$ipni_id), all_h_pow$fq_id),
  function(x) {
    Sys.sleep(.5)
    try(taxize::pow_lookup(x, include = c("distribution")))
  },
  cl = cl, USE.NAMES = FALSE
)

wfo_h_pow_add <- wfo_h_pow_l_add[sapply(wfo_h_pow_l_add, length) == 2] %>%
  lapply(., function(x) x$meta) %>%
  lapply(., rlist::list.flatten) %>%
  lapply(., function(x) x[sapply(x, length) == 1]) %>%
  data.table::rbindlist(., fill = TRUE) %>%
  janitor::clean_names() %>%
  distinct()

wfo_h_pow %<>% bind_rows(wfo_h_pow_add)
write_rds(wfo_h_pow, file = paste0("wfo_h_pow_", lubridate::today(), ".rds"))

# parse out parent information from POWO -----------------------------------------------------------
pow_p <- wfo_h_pow %>%
  mutate(
    parent = case_match(
      hybrid_formula,
      "Hemionitis acrostica H. pteridioides" ~ "H. acrostica × H. pteridioides",
      "S. bulbosum S. tuberosum" ~ "S. bulbosum × S. tuberosum",
      "E. denticulatum and E. myrmecophorum" ~ "E. denticulatum × E. myrmecophorum",
      "S. ferrugineus S. furcatus" ~ "S. ferrugineus × S. furcatus",
      "H. lucidula H. porophila" ~ "Huperzia lucidula × H. porophila",
      "A. chamaejasme × .A obtusifolia" ~ "A. chamaejasme × A. obtusifolia",
      "O;. campestris × O. lapponica" ~ "O. campestris × O. lapponica",
      .default = hybrid_formula
    ) %>%
      str_remove_all(pattern = "\\(|\\)") %>%
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
      str_trim(),
    parent = case_match(
      parent,
      "Erica , pellucida" ~ "Erica pellucida",
      "Euphorbia esula subsp. tomassiniana" ~ "Euphorbia esula subsp. tommasiniana",
      "Smithianthe zebrina" ~ "Smithiantha zebrina",
      "Athyrium foliosum" ~ "Athyrium foliolosum",
      "Athyrium fimnriatum" ~ "Athyrium fimbriatum",
      "Agrimonia pilos" ~ "Agrimonia pilosa",
      "Platanthera biflora" ~ "Platanthera bifolia",
      "Eucalyptus mecrotheca" ~ "Eucalyptus microtheca",
      "Dactyloglossum urvilleana" ~ "Dactylorhiza urvilleana",
      "Cypripedium villosum var. boxallii" ~ "Paphiopedilum villosum var. boxallii",
      "Paphiopedilum villosum var. boxalii" ~ "Paphiopedilum villosum var. boxallii",
      "Populus angustifoila" ~ "Populus angustifolia",
      "Populus monilifera" ~ "Populus deltoides subsp. monilifera",
      "Salix aurtia" ~ "Salix aurita",
      "Salix elaeagnus" ~ "Salix elaeagnos",
      "Salix lapponicum" ~ "Salix lapponum",
      "Salix phylicoflia" ~ "Salix phylicifolia",
      "Populus eugenei" ~ "Populus ×canadensis 'Eugenei'",
      "Sonchustenia bupleuroides" ~ "Sonchus bupleuroides",
      "Sonchustenia leptocephalus" ~ "Sonchus leptocephalus",
      "Baccharis elatoides" ~ "Baccharis elaeoides",
      "Alsophila balanaocarpa" ~ "Alsophila balanocarpa",
      "Polystichum mfibrillosopaleaceum" ~ "Polystichum fibrillosopaleaceum",
      "Platygyria jakonensis" ~ "Lepisorus jakonensis",
      "Asplenophyllitis onopteris" ~ "Asplenium onopteris",
      "Asplenophyllitis scolopendrium" ~ "Asplenium scolopendrium",
      "Moranopteris micropteris" ~ "Moranopteris microlepis",
      "Laeliocattleya purpurata" ~ "Cattleya purpurata",
      "Laeliocattleya tigrina" ~ "Cattleya  tigrina",
      "Anthurium andreanum" ~ "Anthurium andraeanum",
      "Limodorum palustris" ~ "Epipactis palustris",
      "Miltonia warscewiczii var. weltonii" ~ "Miltonia warszewiczii",
      "Oberonia forcipera" ~ "Oberonia forcipifera",
      "Ophrys umblicata subsp. attica" ~ "Ophrys umbilicata subsp. attica",
      "Orchicoeloglossum viride" ~ "Dactylorhiza viride",
      "Orchiserapias mascula" ~ "Orchis mascula",
      "Phlomis gourgaei" ~ "Phlomis bourgaei",
      "Rondeletia miraflorens" ~ "Rondeletia miraflorensis",
      "Scutellaria haematochloa" ~ "Scutellaria chaematochlora",
      "Sideritis ibanyezii" ~ "Sideritis ibanezii",
      "Sophrocattleya crispa" ~ "Cattleya crispa",
      "Sophrocattleya velutina" ~ "Cattleya velutina",
      "Sophrocattleya amethystoglossa" ~ "Cattleya amethystoglossa",
      "Sophrocattleya grandis" ~ "Cattleya  grandis",
      "Teucrium chamsedrys" ~ "Teucrium chamaedrys",
      "Thymus caespititus" ~ "Thymus caespititius",
      "Mahoberberis candidula" ~ "Berberis candidula",
      "Mahoberberis sargentiana" ~ "Berberis sargentiana",
      "Dioscorea cieneyensis" ~ "Dioscorea cienegensis",
      "Procopiphytum davisii subsp. icaricum" ~ "Symphytum davisii subsp. icaricum",
      "Rumex dentatus subsp. klotzschianusrispus" ~ "Rumex dentatus subsp. klotzschianus",
      "Rumex dentatus subsp. kotschyanus" ~ "Rumex dentatus subsp. klotzschianus",
      "Rhododendron smirnovii" ~ "Rhododendron smirnowii",
      "Eleocharis pellucina" ~ "Eleocharis pellucida",
      "Verbascum hassknechtii" ~ "Verbascum haussknechtii",
      "Fimbristylis microcarpa" ~ "Fimbristylis microcarya", # ?
      "Sericobonia ghiesbreghtiana" ~ "Justicia ghiesbreghtiana",
      "Sinocalycalycanthus chinensis" ~ "Calycanthus chinensis",
      "Sinocalycalycanthus floridus" ~ "Calycanthus floridus",
      "Limonium auriculiursifolium" ~ "Limonium auriculae-ursifolium",
      "Macludrania tricuspidata" ~ "Maclura tricuspidata",
      "Myrica caroiniensis" ~ "Myrica caroliniensis",
      "Neoregeila bahiana var. viridis" ~ "Neoregelia bahiana var. viridis",
      "Saxifraga luteovirdis" ~ "Saxifraga luteoviridis",
      "Scabiosa drymeia" ~ "Knautia drymeia",
      "Soldanella hungrica" ~ "Soldanella hungarica",
      "Camellia japnica" ~ "Camellia japonica",
      "Arenaria montana subsp. intricara" ~ "Arenaria montana subsp. intricata",
      "Aconitum bucoviense" ~ "Aconitum bucovinense",
      "Scirpus prolifera" ~ "Scirpus prolifer",
      "Baptisia bracteata var. leucophlaea" ~ "Baptisia bracteata var. leucophaea",
      "Camellia sasanque" ~ "Camellia sasanqua",
      "Camellia saluensis" ~ "Camellia saluenensis",
      "Tradescantia subasper" ~ "Tradescantia subaspera",
      "Centranthus lecoquii" ~ "Centranthus lecoqii",
      "Jamesbrittenia brevifolia" ~ "Jamesbrittenia breviflora",
      "Chaenostoma pristisepala" ~ "Jamesbrittenia pristisepala",
      "Chenopodium berlandieri var. zchackei" ~ "Chenopodium berlandieri var. zschackei",
      "Chenopodium pedunculre" ~ "Chenopodium pedunculare",
      "Cistus lanifer" ~ "Cistus ladanifer",
      "Citroncirus deliciosa" ~ "Citrus deliciosa",
      "Citroncirus maxima" ~ "Citrus maxima",
      "Citroncirus trifoliata" ~ "Citrus trifoliata",
      "Citrus australe" ~ "Citrus australis",
      "Corispermum marschallianum" ~ "Corispermum marschalli",
      "Crassula cultarata" ~ "Crassula cultrata",
      "Cremnophyla linguifolia" ~ "Cremnophila linguifolia",
      "Draba fladnizansis" ~ "Draba fladnizensis",
      "Draba loiseneurii" ~ "Draba loiseleurii",
      "Draba pruniifolia" ~ "Draba brunifolia",
      "Draba zahlbruckeri" ~ "Draba zahlbruckneri",
      "Epilobium hisutum" ~ "Epilobium hirsutum",
      "Erica nivenii" ~ "Erica nivenia",
      "Erica teralix" ~ "Erica tetralix",
      "Eriostemon hispidus" ~ "Eriostemon hispidulus",
      "Eryngium campestris" ~ "Eryngium campestre",
      "Erythraea erythraea" ~ "Centaurium erythraea",
      "Erythraea littorale" ~ "Centaurium littorale",
      "Lilium candicum" ~ "Lilium candidum",
      "Bystropogon origanifolium var. canariae" ~ "Bystropogon origanifolius var. canariae",
      "Greenonium aureum" ~ "Aeonium aureum",
      "Greenonium glutinosum" ~ "Aeonium glutinosum",
      "Heracleum spondylium" ~ "Heracleum sphondylium",
      "Urceocharis moorei" ~ "Urceolina moorei",
      "Urceocharis sanderi" ~ "Urceolina sanderi",
      "Urceocharis urceolata" ~ "Urceolina urceolata",
      "Gastera bicolor" ~ "Gasteria bicolor",
      "Callitriche conocarpa" ~ "Callitriche cophocarpa", # ?
      "Serpias lingua" ~ "Serapias lingua",
      "Catasetum taguariensis" ~ "Catasetum taguariense",
      "Begonia fuchsioidea" ~ "Begonia fuchsioides",
      "Begonia coraliniifolia" ~ "Begonia caroliniifolia",
      "Aegilotrichum monodurum" ~ "Triticum ×monodurum",
      "Agrositanion longifolius" ~ "Elymus longifolius",
      "Psathyrostachys kronenbergii" ~ "Psathyrostachys kronenburgii",
      "Trisetokoeleria spicata" ~ "Koeleria spicata",
      "Trisetokoeleria subalpestre" ~ "Koeleria subalpestre",
      "Cattleya mantiqueira" ~ "Cattleya mantiqueirae",
      "Copernicia brittonanum" ~ "Copernicia brittonorum",
      "Cyperorchis iridioides" ~ "Cymbidium iridioides",
      "Kunzea recurvata" ~ "Kunzea recurva",
      "Ophrys holosericea subsp. candica" ~
        "Ophrys holoserica subsp. candica (E.Nelson ex Soó) Renz & Taubenheim",
      "Ophrys holosericea subsp. grandiflora" ~
        "Ophrys holoserica subsp. grandiflora (H.Fleischm. & Soó) Faurh.",
      "Ophrys holosericea subsp. chestermanii" ~ "Ophrys holoserica subsp. chestermanii J.J.Wood",
      "Ophrys holosericea subsp. oxyrrhynchos" ~
        "Ophrys holoserica subsp. oxyrrhynchos (Tod.) H.Sund.",
      "Anemone vitifolia var. vitifolia" ~ "Anemone vitifolia",
      "Ophrys holosericea subsp. biancae" ~
        "Ophrys holoserica subsp. biancae (Tod.) Faurh. & H.A.Pedersen",
      "Ophrys holosericea subsp. lacaitae" ~ "Ophrys holoserica subsp. lacaitae (Lojac.) W.Rossi",
      "Ophrys holosericea subsp. apulica" ~
        "Ophrys holoserica subsp. apulica (O.Danesch & E.Danesch) Buttler",
      "Ophrys holosericea subsp. parvimaculata" ~
        "Ophrys holoserica subsp. parvimaculata (O.Danesch & E.Danesch) O.Danesch & E.Danesch",
      "Paphiopedilum javanicum var. virens" ~ "Paphiopedilum javanicum",
      "Thymus pulegioides subsp. marschallianus" ~ "Thymus marschallianus",
      "Paphiopedilum philippinense var. roebellinii" ~
        "Paphiopedilum philippinense var. roebelenii (A.H.Kent) P.J.Cribb",
      .default = parent
    )
  ) %>%
  filter(str_count(parent, "\\S+") > 1)

gbif_p <- gbif.cmpfn(
  unique(tmp$parent),
  seq_along(unique(tmp$parent))
) %>%
  filter(str_detect(scientificName, "×| x ") | confidence > 95 | confidence == 0) %>%
  group_by(name) %>%
  filter(n() == 1) %>%
  ungroup()
setdiff(unique(tmp$parent), gbif_p$name)


wfo0_ <- WFO.match(hybr_ecoflora$edited,
                   WFO.data = wfo_data, Fuzzy.one = FALSE, Fuzzy.two = FALSE,
                   verbose = FALSE
)
wfo0_h <- WFO.match(hybr_ecoflora$edited,
                    WFO.data = hybr_wfo, Fuzzy.one = FALSE, Fuzzy.two = FALSE,
                    verbose = FALSE
)
wfo_ <- wfo0_ %>%
  filter(taxonID != "", Hybrid == "×") %>%
  select(spec.name.ORIG, taxonID, scientificName, source) %>%
  distinct() %>%
  group_by(spec.name.ORIG) %>%
  filter(n() == 1)

# get POW taxon names and IDs based on canonical name matches
out1_ <- pbapply::pbsapply(
  unique(hybr_ecoflora$edited3),
  function(x) {
    Sys.sleep(.1)
    try(taxize::get_pow_(x, accepted = TRUE, messages = FALSE))
  },
  cl = cl, USE.NAMES = FALSE
)

pow_ <- taxize::get_pow(unique(hybr_ecoflora$edited3))

tmp <- all_h %>%
  filter(is.na(ipni_id)) %>%
  pull(name_edited)

gbif_p <- gbif.cmpfn(unique(tmp), seq_along(unique(tmp))) %>%
  filter(str_detect(scientificName, "×| x ") | confidence > 95 | confidence == 0) %>%
  group_by(name) %>%
  filter(n() == 1) %>%
  ungroup()
setdiff(tmp, gbif_p$name)

wfo_h$name_ipni
tmp <- wfo_h %>% filter(scientific_name_id %in% setdiff(scientific_name_id, ipni_h$id))


# ------
# check if all IPNI IDs provided in WFO are valid --------------------------------------------------
wfo_h0_id <- unique(wfo_h0$ipni_id[!is.na(wfo_h0$ipni_id)])
pow_info_l <- pbapply::pbsapply(
  wfo_h0_id,
  function(x) {
    Sys.sleep(.2)
    try(taxize::pow_lookup(x))
  },
  cl = cl, USE.NAMES = FALSE
)

wfo_h0_pow <- pow_info_l[sapply(pow_info_l, length) == 2] %>%
  lapply(., function(x) x$meta) %>%
  lapply(., rlist::list.flatten) %>%
  lapply(., function(x) x[sapply(x, length) == 1]) %>%
  data.table::rbindlist(., fill = TRUE) %>%
  janitor::clean_names() %>%
  distinct()

pow_info_l_add <- pbapply::pbsapply(
  setdiff(wfo_h0_id, wfo_h0_pow$fq_id),
  function(x) {
    Sys.sleep(.5)
    try(taxize::pow_lookup(x))
  },
  cl = cl, USE.NAMES = FALSE
)

wfo_h0_pow_add <- pow_info_l_add[sapply(pow_info_l_add, length) == 2] %>%
  lapply(., function(x) x$meta) %>%
  lapply(., rlist::list.flatten) %>%
  lapply(., function(x) x[sapply(x, length) == 1]) %>%
  data.table::rbindlist(., fill = TRUE) %>%
  janitor::clean_names() %>%
  distinct()

wfo_h0_pow %<>% bind_rows(wfo_h0_pow_add)
rm(wfo_h0_pow_add)
setdiff(wfo_h0_id, wfo_h0_pow$fq_id)