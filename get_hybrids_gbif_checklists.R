# Marina Golivets
# BiCIKL Hackathon (Sept 20-24, 2021, Meise)
# Extract hybrid and parent taxon names from checklists that are available via GBIF


# load packages ------------------------------------------------------------------------------

pkgs <- c("here", "tidyverse", "magrittr")
sapply(pkgs, require, character.only = TRUE)


# initialize parallelization ------------------------------------------------------------------

no_cores <- parallel::detectCores()
cl <- parallel::makeCluster(no_cores)


# get checklists using GBIF API ---------------------------------------------------------------

# code from: 
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

write_csv2(out, here("data/checklists_gbif_raw.csv"))


# get hybrid names from checklists in GBIF ------------------------------------------------------

gbif_chl <- data.table::fread(here("data/checklists_gbif_raw.csv")) %>%
  filter(kingdom == "Plantae" | is.na(kingdom)) %>%
  pivot_longer(c(scientificName, accepted), values_to = "verbatim") %>%
  filter(!is.na(verbatim)) %>%
  select(kingdom, verbatim) %>%
  filter(!grepl(" '", verbatim)) %>% # remove cultivars
  filter(!grepl("phytoplasma|mycoplasma|Bacillus|Mörkt|Fläcknycklar|Asterlateriflorus",
    verbatim, ignore.case = TRUE
  )) %>%
  filter(!grepl(
    "^[a-z]\\. |^[a-z]{2}\\. |symbiont|Betacoccus|Yucatan|Drosophila|Mandarine|Autumn Joy",
    verbatim, ignore.case = TRUE
  )) %>%
  filter(!(verbatim %in% c(
    "Corallinales sp. X_Twist-2019", "×Hulthemosa Juz. (Hulthemia × Rosa)",
    "Maranta pulchella (Jacob-Makoy x É.Morren) É.Morren",
    "Campsis radicans (L.) Seem. (incl. C. x tagliabuana (Vis.) Rehd.)",
    "Pseudochantransia serpens (Israelson) Sheathe x Kumano",
    "Mouriria doriana Saldanha x Cogn.", "Hyptis gratiolaefolia A. St-Hil. x Benth.",
    "Calyplectus adenophyllus Schott x Mart.", "Apicomplexa sp. X STN-2017",
    "×Conyzigeron Rauschert", "×Triticosecale Wittm. ex A.Camus",
    "×Dactylodenia Garay & H.R.Sweet", "×Dactyloglossum P.F.Hunt & Summerh.",
    "×Gymnanacamptis Asch. & Graebn.", "×Orchiaceras E.G.Camus", "×Gymnaglossum Rolfe",
    "×Pseudadenia P.F.Hunt", "×Pseudorhiza P.F.Hunt", "×Tripleurothemis Stace",
    "×Crataemespilus E.G.Camus", "×Heucherella H.R.Wehrh.", "×Cupressocyparis Dallim.",
    "×Agropogon P.Fourn.", "×Calammophila Brand", "×Elytrordeum Hyl.", 
    "×Festulolium Asch. & Graebn.", "×Festulpia Melderis ex Stace & R.Cotton",
    "× Gladanthera J.M.Wright × Homoglad Ingram")
  )) %>% # remove taxa which are not species or hybrids
  mutate(edited = gsub('"', '', verbatim)) %>%
  mutate(edited = gsub("\\(\\?\\)|\\?", "", edited)) %>%
  mutate(edited = gsub(" x |^x |×| × ×", " × ", edited, ignore.case = TRUE)) %>%
  mutate(edited = gsub("\\(x ", "(× ", edited, ignore.case = TRUE)) %>%
  mutate(edited = gsub("_x ", " ", edited, ignore.case = TRUE)) %>%
  mutate(edited = str_trim(str_squish(edited))) %>%
  mutate(edited = gsub("Asterales ", "Aster ", edited)) %>%
  rowwise() %>%
  mutate(genus = str_split(gsub("^× ", "", edited), pattern = " ", simplify = TRUE)[1]) %>%
  mutate(genus = str_to_sentence(genus)) %>%
  distinct()


# match genus names to the GBIF backbone
genera <- gbif_chl %>%
  pull(genus) %>%
  unique() 
genera %<>%  pbapply::pblapply(rgbif::name_backbone, rank = "genus", cl = cl) %>%
  bind_rows() %>%
  mutate(genus = genera) %>%
  filter(kingdom == "Plantae" | is.na(kingdom))

# clean data further
hybrids_gbif_chl <- gbif_chl %>%
  filter(genus %in% genera$genus) %>% # remove non-plant entries
  select(verbatim, edited) %>%
  mutate(edited = gsub(
    " NULL|\\(f × m\\)|\\(female × male\\)", "", edited)) %>%
  mutate(edited = gsub(
    "pro\\.sp\\.|\\(pro sp\\.)|, perhaps.*$|L\\. × L\\. × L\\.", "", edited)) %>%
  mutate(edited = gsub("ssp\\. ", "subsp. ", edited)) %>%
  mutate(edited = gsub(" forma ", " f. ", edited)) %>%
  mutate(edited = gsub(" sens.lat.", " s. l.", edited)) %>%
  mutate(edited = gsub(" var\\. ×", " var.", edited)) %>%
  mutate(edited = gsub(" var ", " var. ", edited)) %>%
  mutate(edited = gsub(" f ", " f. ", edited)) %>%
  mutate(edited = gsub(" subsp\\. ×", " subsp.", edited)) %>%
  mutate(edited = gsub(" f\\. ×", " f.", edited)) %>%
  mutate(edited = gsub(" physcomitrium ", " Physcomitrium ", edited)) %>%
  mutate(edited = gsub(" pohlia ", " Pohlia ", edited)) %>%
  mutate(edited = gsub(" ditrichum ", " Ditrichum ", edited)) %>%
  mutate(edited = gsub(" pleuridium ", " Pleuridium ", edited)) %>%
  mutate(edited = gsub(" funaria ", " Funaria ", edited)) %>%
  mutate(edited = gsub(" voitia ", " Voitia ", edited)) %>%
  mutate(edited = gsub("Pilosell ", "Pilosella ", edited)) %>%
  mutate(edited = gsub("^salix ", "Salix ", edited)) %>%
  mutate(edited = gsub("filix femina ", "filix-femina ", edited)) %>%
  mutate(edited = gsub(" filix mas| felix-mas", " filix-mas", edited)) %>%
  mutate(edited = gsub(" ruta muraria", " ruta-muraria", edited)) %>%
  mutate(edited = gsub("novi‚Äìbelgii", "novi-belgii", edited)) %>%
  mutate(edited = gsub(" cystopteris ", " Cystopteris ", edited)) %>%
  mutate(edited = gsub(" asplenium ", " Asplenium ", edited)) %>%
  mutate(edited = gsub(" stipa ", " Stipa ", edited)) %>%
  mutate(edited = gsub(" dryopteris ", " Dryopteris ", edited)) %>%
  mutate(edited = gsub(" adiantum ", " Adiantum ", edited)) %>%
  mutate(edited = gsub("^Festulpia Festuca |^x Festulolium Festuca ", "Festuca ", edited)) %>%
  mutate(edited = gsub(" stoloniferar", " stolonifera", edited)) %>%
  mutate(edited = gsub(" Crocosmiiflora", " crocosmiiflora", edited)) %>%
  mutate(edited = gsub(" KotukhovS. orientalis ", " Kotukhov = S. orientalis ", edited)) %>%
  mutate(edited = gsub("Festuca × Festulolium", "Festulolium", edited)) %>%
  mutate(edited = case_when(
    verbatim == "Rorippa podolica x Zapal." ~ "Rorippa × podolica Zapal.",
    verbatim == "Rorippa viaria x Zapal." ~ "Rorippa × viaria Zapal.",
    verbatim == "Rorippa sodalis x Zapal." ~ "Rorippa × sodalis Zapal.",
    verbatim == "Rorippa stenophylla x Borbás ex Nyman" ~ "Rorippa × stenophylla Borbás ex Nyman",
    verbatim == "Ranunculus x bachii sens.lat. (R. aquatilis or trichophyllus x fluitans)" ~
    "Ranunculus x bachii s. l.",
    verbatim == "Filipendula camtschatica x sp. (F. x purpurea)" ~
    "Filipendula × camtschatica (F. × purpurea)",
    verbatim == "Draba nichanaica x O.E.Schulz" ~ "Draba × nichanaica O.E.Schulz",
    verbatim == "Draba pseudonivalis x N.Busch" ~ "Draba × pseudonivalis N.Busch",
    verbatim == "Spiraea salicifolia (x billardii)" ~ "Spiraea salicifolia × billardii",
    verbatim == "Potamogeton x angustifolius (x zizii)" ~ "Potamogeton × angustifolius × zizii",
    TRUE ~ edited
  )) %>%
  mutate(edited = gsub(
    "Erigeron × Conyza |Elymus × Hordeum |Peltaria × Bornmuellera ", "", edited
  )) %>%
  mutate(edited = gsub("Dactylorhiza × Gymnadenia |Leymus × Elymus|\\/Xiphinema.*$", "", edited
  )) %>%
  mutate(edited = gsub(
    "Aceras × Orchis |Coeloglossum × Dactylorhiza |Peltaria × Bornmuellera ", "", edited
  )) %>%
  filter(!grepl("^[A-Z][a-z]+ × [A-Z][a-z][a-z]+", edited)) %>%
  mutate(edited = str_trim(str_squish(edited))) %>%
  mutate(edited = gsub(" or ", "=", edited)) %>%
  rowwise() %>%
  mutate(a = str_extract(edited, "\\(.* × .*\\)|= .* × .*|= × .*|\\(× .*\\)|\\[.* × .*\\]")) %>%
  mutate(a = gsub("\\(|\\)|=|\\[|\\]", "", a)) %>%
  mutate(a = gsub("^.*Camus", "", a)) %>%
  mutate(edited = gsub("\\(.* × .*\\)|= .* × .*|= × .*|\\(× .*\\)|\\[.* × .*\\]", "", edited)) %>%
  mutate(across(edited:a, ~ str_trim(str_squish(.)))) %>%
  # edit the names that were messed up because of nested parentheses 
  mutate(edited = case_when(
    verbatim == "Rumex aquaticus x (aquaticus x hydrolapathum)" ~
    "Rumex aquaticus × (aquaticus × hydrolapathum)",
    verbatim == "Mentha aquatica x (aquatica x spicata)" ~
    "Mentha aquatica × (aquatica × spicata)",
    verbatim == "Dactylorhiza incarnata (L.) Soó × Dactylorhiza maculata (L.) Soó" ~
    "Dactylorhiza incarnata (L.) Soó × Dactylorhiza maculata (L.) Soó",
    edited == "× Festulolium braunii" ~ "× Festulolium braunii (K. Richt.) A. Camus",
    verbatim == "Spiranthes cernua (L.) Richard x S. odorata  (Nuttall) Lindley" ~
    "Spiranthes cernua (L.) Richard × S. odorata  (Nuttall) Lindley",
    TRUE ~ edited
  )) %>%
  mutate(a = ifelse(
    edited %in% c(
      "Rumex aquaticus × (aquaticus × hydrolapathum)",
      "Mentha aquatica × (aquatica × spicata)",
      "Dactylorhiza incarnata (L.) Soó × Dactylorhiza maculata (L.) Soó",
      "Spiranthes cernua (L.) Richard × S. odorata (Nuttall) Lindley"
    ),
    NA, a
  )) %>%
  mutate(genus = str_split(gsub("^× ", "", edited), pattern = " ", simplify = TRUE)[1]) %>%
  mutate(a = ifelse(grepl("^[a-z]", a), paste(genus, a), a)) %>%
  mutate(a = gsub("^[A-Za-z]\\.|^[A-Za-z] ", paste(genus, " "), a)) %>%
  mutate(across(edited:a, ~ str_trim(str_squish(.)))) %>%
  pivot_longer(edited:a, values_to = "edited") %>%
  mutate(edited = gsub("^\\(| \\)$|,$|\\(\\)$", "", edited)) %>%
  mutate(across(edited, ~ str_trim(str_squish(.)))) %>%
  select(-genus, -name) %>%
  filter(!is.na(edited)) %>%
  mutate(parent = str_split(edited, pattern = "×")) %>%
  unnest("parent") %>%
  mutate(parent = str_trim(str_squish(parent))) %>%
  filter(parent != "") %>%
  rowwise() %>%
  mutate(genus = str_split(
    gsub("^× ", "", edited),
    pattern = " ", simplify = TRUE
  )[1]) %>%
  mutate(genus = str_to_sentence(genus)) %>%
  mutate(across(parent:genus, ~ str_trim(str_squish(.)))) %>%
  filter(parent != genus) %>%
  mutate(parent = gsub("^\\(|\\)$|,$|\\(\\)$", "", parent)) %>%
  mutate(across(parent, ~ str_trim(str_squish(.)))) %>%
  mutate(parent = gsub("^[A-Za-z]\\.|^[A-Za-z] ", paste(genus, " "), parent)) %>%
  mutate(parent = ifelse(grepl("^[a-z]", parent) & !grepl("^subsp\\.", parent),
    paste(genus, parent), parent
  )) %>%
  mutate(parent = ifelse(
    grepl("^subsp\\.", parent),
    paste(
      str_split(gsub("^× ", "", edited), pattern = " ", simplify = TRUE)[1],
      str_split(gsub("^× ", "", edited), pattern = " ", simplify = TRUE)[2], parent
    ), parent
  )) %>%
  mutate(across(parent, ~ str_trim(str_squish(.)))) %>%
  mutate(parent = ifelse(parent == edited, NA, parent)) %>%
  mutate(parent = ifelse(gsub(" × ", " ", edited) == parent, NA, parent)) %>%
  mutate(across(parent, ~ str_trim(str_squish(.)))) %>%
  select(-genus) %>%
  distinct()

hybrids_gbif_chl %<>%
  select(verbatim, edited) %>%
  full_join(
    hybrids_gbif_chl %>%
    select(verbatim, parent) %>%
      filter(!is.na(parent)) 
    ) %>%
  distinct() %>%
  write_csv2(here("data/hybrids_gbif_checklists.csv"))


# source(here("fn_match_names_to_gbif.R"))
# gbif.cmpfn <- compiler::cmpfun(gbif.fn)
# parent_matched <- gbif.cmpfn(
#   taxonName = unique(hybrids_gbif_chl$parent),
#   taxonID = seq_along(unique(hybrids_gbif_chl$parent))
# ) %>%
#   rename(parent = taxonName)
# hybrids_gbif_chl %<>%
#   left_join(parent_matched)
# 
# 
# hybrids_gbif_chl %>%
#   filter(is.na(kingdom)) %>%
#   pull(parent) %>%
#   unique()

parallel::stopCluster(cl)