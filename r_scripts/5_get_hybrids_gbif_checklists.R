# Marina Golivets
# BiCIKL Hackathon (Sept 20-24, 2021, Meise, Belgium)
# Extract hybrid and parent taxon names from checklists that are available via GBIF


# get checklists using GBIF API --------------------------------------------------------------------
# R code from:
# https://github.com/PietrH/1-lost_parents/blob/main/src/getting_hybrids_from_checklists_gbif.R

limit <- 1000
offset <- 0
out <- tibble()
for (offset in seq(0, 40000, by = limit)) {
  print(offset)

  url <-
    sprintf(
      "https://api.gbif.org/v1/species/search?name_type=HYBRID&limit=%i&offset=%i",
      limit,
      offset
    )
  response <- httr::GET(url)
  parsed <- httr::content(response, as = "parsed")
  parsed_columns <-
    lapply(parsed$results, function(x) {
      x[c(
        "key",
        "accepted",
        "kingdom",
        "scientificName",
        "speciesKey",
        "taxonID"
      )]
    })
  out <- bind_rows(out, data.table::rbindlist(parsed_columns, fill = TRUE)) %>%
    janitor::remove_empty(which = "cols")
  message(sprintf("we have %i rows in our out with offset %i", nrow(out), offset))
}

data.table::setDT(out)
data.table::setkey(out, speciesKey)

write_csv2(out, here(paste0("output/checklists_gbif_raw_", lubridate::today(), ".csv")))

# get hybrid names from checklists in GBIF ---------------------------------------------------------

gbif_chl0 <- data.table::fread(
  here("output/checklists_gbif_raw_2023-05-12.csv")
) %>%
  filter(kingdom == "Plantae" | is.na(kingdom)) %>% # filter out non-plants
  select(scientificName, accepted) %>%
  mutate(
    accepted = if_else(
      str_replace_all(scientificName, " x ", " × ") == accepted,
      NA_character_, accepted, missing = accepted
    )
  ) %>%
  distinct() %>%
  rowid_to_column("id") %>%
  group_by(scientificName)

# leave only matches with a single accepted name
gbif_chl1 <- bind_rows(
  gbif_chl0 %>%
    filter(n() == 1),
  gbif_chl0 %>%
    filter(n() > 1) %>%
    filter(!is.na(accepted))
) %>%
  ungroup() %>%
  glimpse()

# str_to_exclude <- c(
#   "'", "Autumn Joy", # cultivars
#   "incl\\. ", # more than one taxon
#   "^[A-Z][a-z]+ x [A-Z][a-z.]+$", "^[A-Z][a-z]+ × [A-Z][a-z.]+$",
#   "^×[A-Z][a-z]+ [A-Z]", "^[A-Za-z]+ [A-Za-z.]+$", # genera
#   "plasma", "kungsljus", "knycklar", "Mandarine", "Bacillus",
#   "symbiont", "Drosophila", "X_Twist-2019", "X STN-2017",
#   "Jacob-Makoy x É.Morren"
# )
#
# gbif_chl %<>%
#   filter(
#     str_detect(scientificName,
#       str_c(str_to_exclude, collapse = "|"),
#       negate = TRUE
#     )
#   ) %>%
#   mutate(accepted = case_when(
#     str_detect(accepted, str_c(str_to_exclude, collapse = "|")) ~ NA_character_,
#     TRUE ~ accepted
#   ))

gbif_chl2 <- gbif_chl1 %>%
  pivot_longer(c(scientificName, accepted), values_to = "verbatim", values_drop_na = TRUE) %>%
  select(-name) %>%
  distinct() %>%
  mutate(
    verbatim,
    edited = str_remove_all(
      verbatim,
      str_c(c('"', "\\[.*\\]", "\\(\\?\\)|\\?", "_x |_X"," -$"), collapse = "|")
    ) %>%
      str_replace_all(" (?i)x |^(?i)x |×| × ×", " × ") %>%
      str_replace_all("\\((?i)x ", "(× ") %>%
      str_squish() %>%
      str_trim() %>%
      str_replace_all("Asterales ", "Aster "),
    genus = str_remove(edited, "^× ") %>%
      str_extract("\\w*") %>%
      str_to_sentence()
  ) %>%
  filter(
    str_detect(
      edited, 
      "^×[A-Z][a-z]+ [A-Z]|^[A-Za-z]+ × [A-Za-z]+$|^[A-Z][a-z]+ [A-Z]|incl\\.|fläcknycklar",
      negate = TRUE
    )
  ) %>%
  distinct()

# match genus names to the GBIF backbone
genera0 <- gbif_chl2 %>%
  pull(genus) %>%
  unique()
genera <- genera0 %>% 
  pbapply::pblapply(rgbif::name_backbone, rank = "genus", cl = cl) %>%
  bind_rows() %>%
  filter(
    kingdom == "Plantae" | is.na(kingdom) | verbatim_name == "Ammophila",
    !(verbatim_name %in% c("Green", "Bacillus", "Mörkt", "Yucatan"))
  )
setdiff(genera0, genera$verbatim_name)
gbif_chl3 <- gbif_chl2 %>%
  filter(genus %in% genera$verbatim_name) %>% # remove non-plants
  select(-genus)

# clean data further
gbif_chl4 <- gbif_chl3 %>%
  mutate(
    edited = edited %>%
      str_remove_all(pattern = str_c(c(" NULL", "[0-9]"), collapse = "|")) %>%
      str_replace_all(
        c(
          "^× " = "",
          "ssp\\. " = "subsp. ",
          " forma " = " f. ",
          " sens.lat." = " s. l.",
          " var\\. ×" = " var.",
          " var " = " var. ",
          " f " = " f. ",
          " subsp\\. ×" = " subsp.",
          " f\\. ×" = " f.",
          " physcomitrium " = " Physcomitrium ",
          " pohlia " = " Pohlia ",
          " physcomitrium " = " Physcomitrium ",
          " pohlia " = " Pohlia ",
          " ditrichum " = " Ditrichum ",
          " pleuridium " = " Pleuridium ",
          " funaria " = " Funaria ",
          " voitia " = " Voitia ",
          "Pilosell " = "Pilosella ",
          "^salix " = "Salix ",
          "filix femina " = "filix-femina ",
          " filix mas| felix-mas" = " filix-mas",
          " ruta muraria" = " ruta-muraria",
          "novi‚Äìbelgii" = "novi-belgii",
          " cystopteris " = " Cystopteris ",
          " asplenium " = " Asplenium ",
          " stipa " = " Stipa ",
          " dryopteris " = " Dryopteris ",
          " adiantum " = " Adiantum ",
          " stoloniferar" = " stolonifera",
          " Crocosmiiflora" = " crocosmiiflora",
          "Festuca × Festulolium" = "Festulolium",
          "Festulolium Festuca" = "Festuca"
        )
      )
  ) %>%
  mutate(
    edited = case_when(
      verbatim == "Spiraea salicifolia (x billardii)" ~ "Spiraea salicifolia × billardii",
      verbatim == "Potamogeton x angustifolius (x zizii)" ~ "Potamogeton × angustifolius × zizii",
      verbatim == "Dactylorhiza incarnata (L.) Soó × Dactylorhiza maculata (L.) Soó" ~
        "Dactylorhiza incarnata (L.) Soó × Dactylorhiza maculata (L.) Soó",
      str_detect(edited, "Festulolium braunii") ~ "Festulolium braunii (K. Richt.) A. Camus",
      str_detect(edited, "Maranta pulchella ") ~ NA_character_,
      verbatim == "Spiranthes cernua (L.) Richard x S. odorata  (Nuttall) Lindley" ~
        "Spiranthes cernua × S. odorata",
      str_detect(edited, "Alcea ficifolia") ~ "Alcea rosea × rugosa",
      str_detect(edited, "Stipa ladakhensis") ~ "Stipa klimesii × S. purpurea s. l.",
      str_detect(edited, "Mentha × piperita L. nsubsp. piperita") ~
        "Mentha aquatica L. × spicata subsp. glabrata",
      str_detect(edited, "Trichophorum cespitosum subsp. nothofoersteri") ~
        "Trichophorum cespitosum × germanicum",
      str_detect(edited, "Sorbus mougeotii") ~ "Sorbus aucuparia × aria",
      str_detect(edited, "Chara papillosa") ~ "Chara hispida × contraria",
      TRUE ~ edited
    ) %>%
      str_remove_all(
        str_c(
          c(
            "Erigeron × Conyza ", "Elymus × Hordeum ", "Peltaria × Bornmuellera ",
            "Dactylorhiza × Gymnadenia ", "Leymus × Elymus", "\\/Xiphinema.*$",
            "Aceras × Orchis ", "Coeloglossum × Dactylorhiza ", "Peltaria × Bornmuellera "
          ),
          collapse = "|"
        )
      )
  ) %>%
  mutate(type = case_when(
    str_count(edited, "×") > 1 ~ "formula",
    str_detect(edited, "^[A-Za-z-]+ ×|^×|!×") ~ "name",
    TRUE ~ "formula"
  ))

gbif_chl5 <- gbif_chl4 %>%
  filter(type == "formula") %>%
  transmute(
    formula = edited,
    parent = str_split(formula, pattern = " × ")
  ) %>%
  unnest(cols = "parent") %>% # parse out hybrid parents from hybrid formulas
  mutate(across(everything(), ~ str_squish(.) %>% str_trim())) %>%
  distinct() %>%
  mutate(
    genus = str_extract(formula, "\\w*"),
    parent = str_replace(parent, "^[A-Z]\\.", genus) %>%
      if_else(str_detect(., "^subsp |^subsp\\."), paste(str_extract(formula, "\\w* \\w*"), .), .) %>%
      if_else(str_detect(., "^[a-z]"), paste(genus, .), .) %>%
      if_else(str_detect(., " ", negate = TRUE), NA_character_, .)
  ) %>%
  filter(!is.na(parent)) %>%
  mutate(parent = str_remove(parent, "\\(|\\)")) %>%
  group_by(formula) %>%
  filter(n() == 2)

gbif_chl6 <- gbif_chl5 %>%
  left_join(gbif_chl4 %>%
    select(id, edited, verbatim) %>%
    distinct(),
  by = c("formula" = "edited")
  ) %>%
  left_join(
    gbif_chl4 %>%
      filter(type == "name") %>%
      transmute(name = edited, id)
  ) %>%
  select(-id, -genus) %>%
  mutate(
    name = case_when(
      str_detect(formula, "Alcea rosea × rugosa") ~ "Alcea ficifolia",
      str_detect(formula, "Stipa klimesii × S. purpurea s. l.") ~ "Stipa ladakhensis",
      str_detect(formula, "Mentha aquatica L. × spicata subsp. glabrata") ~
        "Mentha × piperita L. subsp. piperita",
      str_detect(formula, "Trichophorum cespitosum × germanicum") ~
        "Trichophorum cespitosum subsp. nothofoersteri",
    str_detect(formula, "Sorbus aucuparia × aria") ~ "Sorbus mougeotii",
    str_detect(formula, "Chara hispida × contraria") ~ "Chara papillosa",
    TRUE ~ name
  )) %>%
  pivot_longer(c(name, formula), values_to = "edited", names_to = "type", values_drop_na = TRUE) %>%
  distinct() %>%
  mutate(source = "gbif checklists")

names <- gbif_chl6 %>% 
  filter(type == "name")
formulas <- gbif_chl6 %>% 
  filter(type == "formula", !(verbatim %in% names$verbatim))
gbif_chl7 <- bind_rows(names, formulas) %>%
  mutate(edited = if_else(type == "name", str_replace(edited, " × ", " ×"), edited)) %>%
  write_csv2(here("output/hybrids_gbif_checklists.csv"))
