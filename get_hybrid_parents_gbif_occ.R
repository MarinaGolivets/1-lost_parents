# Marina Golivets
# Clean hybrid names originating from GBIF occurrence data
# Extract parent taxa
# Match taxon names to the GBIF taxonomy


# The list was split into three parts:
# 1) list of hybrid names + hybrid formulas (with "=")
# 2) list of hybrid names
# 3) the rest


# load packages ------------------------------------------------------------------------------------

pkgs <- c("here", "tidyverse", "parallel", "pbapply", "magrittr")
sapply(pkgs, require, character.only = TRUE)


# part 1 -------------------------------------------------------------------------------------------

hybrids_gbif_occ1 <- read_table(
  here("data/hybridnamesWithHybridFormulasFromOccurences.txt"),
  col_names = FALSE,
) %>%
  rename(verbatim = X1) %>%
  mutate(edited = gsub("\\?|- Festul asc", "", verbatim)) %>%
  mutate(edited = gsub(" ssp\\. ", " subsp. ", edited)) %>%
  mutate(edited = gsub("\\.x |\\.x\\.", ". × ", edited)) %>%
  mutate(edited = gsub("×| x | X |^X ", " × ", edited)) %>%
  mutate(across(edited, ~ str_trim(str_squish(.)))) %>%
  rowwise() %>%
  mutate(genus = str_split(edited, pattern = " ", simplify = TRUE)[1]) %>%
  separate(edited, into = c("a", "b"), sep = "=|<=>", remove = FALSE) %>%
  mutate(across(a:b, ~ str_trim(str_squish(.)))) %>%
  mutate(a = gsub("\\($|\\[$", "", a)) %>%
  mutate(b = gsub("^\\(|\\)$|\\]$", "", b)) %>%
  mutate(across(a:b, ~ str_trim(str_squish(.)))) %>%
  mutate(b = ifelse(
    grepl("^[a-zA-Z]\\.|^× [a-zA-Z]\\.|^× [a-z]+ |^[A-Z] |^[a-zA-Z]{2}\\.", b),
    paste(genus, gsub("^[a-zA-Z]\\.|^× [a-zA-Z]\\.|^× [a-z]+ |^[A-Z] |^[a-zA-Z]{2}\\.", "", b)), b
  )) %>%
  rowwise() %>%
  mutate(b = ifelse(
    grepl("× [a-zA-Z]\\.", b), gsub("× [a-zA-Z]\\.", paste("× ", genus), b), b
  )) %>%
  mutate(b = ifelse(grepl("^[a-z]", b), paste(genus, b), b)) %>%
  mutate(a = ifelse
  (grepl("× [a-zA-Z]\\.", a), gsub("× [a-zA-Z]\\.", paste("× ", genus), a), a)) %>%
  mutate(b = gsub("^× ", "", b)) %>%
  mutate(across(a:b, ~ str_trim(str_squish(.)))) %>%
  rowwise() %>%
  mutate(a_f = str_detect(a, paste(genus, "×"))) %>%
  mutate(b_f = str_detect(b, paste(genus, "×"))) %>%
  mutate(c_f = str_detect(a, "×")) %>%
  mutate(d_f = str_detect(b, "×")) %>%
  mutate(formula = case_when(
    a_f == FALSE & b_f == TRUE ~ a,
    b_f == FALSE & a_f == TRUE ~ b,
    b_f == FALSE & a_f == FALSE & c_f == TRUE ~ a,
    b_f == FALSE & a_f == FALSE & d_f == TRUE ~ b
  )) %>%
  mutate(edited = case_when(
    a_f == TRUE & b_f == FALSE ~ a,
    b_f == TRUE & a_f == FALSE ~ b,
    b_f == FALSE & a_f == FALSE & c_f == FALSE ~ a,
    b_f == FALSE & a_f == FALSE & d_f == FALSE ~ b
  )) %>%
  mutate(edited = ifelse(is.na(edited) & !grepl("×", b), b, edited)) %>%
  mutate(formula = ifelse(is.na(formula) & grepl("×", a), a, formula)) %>%
  mutate(edited = ifelse(is.na(edited) & !grepl("×", a), a, edited)) %>%
  mutate(edited = ifelse(formula == "Platanus × hybrida", "Platanus × hybrida", edited)) %>%
  mutate(formula = ifelse(formula == "Platanus × hybrida", NA, formula)) %>%
  mutate(parent = str_split(formula, pattern = "×")) %>%
  unnest("parent") %>%
  mutate(across(parent, ~ str_trim(str_squish(.)))) %>%
  rowwise() %>%
  mutate(genus = str_split(formula, pattern = " ", simplify = TRUE)[1]) %>%
  mutate(parent = case_when(
    grepl("^[a-z]", parent) ~ paste(genus, parent),
    grepl("^[A-Z]\\.", parent) ~
    paste(genus, gsub("^[A-Z]\\.", "", parent)),
    TRUE ~ parent
  )) %>%
  select(verbatim, edited, parent) %>%
  ungroup()


# part 2 -------------------------------------------------------------------------------------------

hybrids_gbif_occ2 <- read_table(here("data/hybridnamesFromOccurences.txt"), col_names = FALSE) %>%
  rename(verbatim = X1) %>%
  filter(!grepl("'", verbatim)) %>%
  filter(!(verbatim %in% c(
    "S.ilicifolius ×S.soldanelloides", "S.multisetus ×S.triangulatus", "pelargonium ×l.h.bailey",
    "L.pilosa ×L.cuneata", "Rumex ×longif.-pseudonatr.", "Abelai ×grandiflora"
  ))) %>%
  mutate(edited = gsub("×", " × ", verbatim)) %>%
  mutate(edited = str_trim(str_squish(edited))) %>%
  mutate(edited = case_when(
    verbatim == "Tripsacum ×maizar-zopilotense-hybrid)" ~ "Tripsacum maizar × zopilotense",
    TRUE ~ edited
  )) %>%
  filter(!str_detect(edited, " × [A-Z][a-z]")) %>%
  filter(!grepl("× sp\\.$|× hybr\\.$|× L\\.$", edited))


# part 3 -------------------------------------------------------------------------------------------

hybrids_gbif_occ3 <- read_csv(here("data/2021-09-23_14-47_hybridnamesUniq2_parents.csv")) %>%
  select(input_strings) %>%
  rename(verbatim = input_strings) %>%
  mutate(edited = str_trim(str_squish(verbatim))) %>%
  filter(
    !grepl(" '|' |´|‘|’| cv\\.| f\\. [A-Z]| var\\. [A-Z]", edited)
  ) %>% # exclude cultivar names
  filter(!grepl("=", edited)) %>% # exclude names with '=' (already processed)
  mutate(edited = gsub("× ×|× x | x x ", "×", edited, ignore.case = TRUE)) %>%
  mutate(edited = gsub(" x |^x |×|x_|^x-", " × ", edited, ignore.case = TRUE)) %>%
  mutate(edited = gsub("[0-9]", "", edited)) %>% # remove numbers
  mutate(edited = str_trim(str_squish(edited))) %>%
  mutate(
    edited = gsub("^[a-z]+ceae ", "", edited, ignore.case = TRUE)
  ) %>% # remove family names
  mutate(edited = gsub(
    "Complex hybrid|doubled|sterile triploid|possible hybrid",
    "", edited
  )) %>%
  mutate(edited = gsub(
    "autotetraploid|diploid|\\(allo|suggested hybrid|form A|form B|unranked",
    "", edited
  )) %>%
  mutate(edited = gsub(".*or possibly", "", edited)) %>%
  mutate(edited = gsub("Cliffortia sp\\.", "Cliffortia", edited)) %>%
  filter(!grepl(
    "× \\?$| × [A-Z][a-z]+ sp\\.|× unknown$| unknown species|indet\\. sp\\.",
    edited,
    ignore.case = TRUE
  )) %>% # exclude hybrid formulas with partially unknown parent
  filter(!grepl(
    "<ALL>|× [a-zA-Z]$|× [a-zA-Z]\\.$|× ssp\\.$| or ",
    edited,
    ignore.case = TRUE
  )) %>% # exclude hybrid formulas with partially unknown parent
  mutate(edited = gsub("unknown", "", edited)) %>%
  mutate(edited = gsub("\\?|_|\\[|\\]| -|- |\\+|-$|--$", "", edited)) %>%
  mutate(edited = gsub("pro sp\\.|pro\\. sp\\.", "", edited, ignore.case = TRUE)) %>%
  filter(!grepl(" sp\\.", edited)) %>%
  mutate(edited = gsub(" ssp\\. ", " subsp. ", edited)) %>%
  mutate(edited = gsub(" nssp\\. ", " nothosubsp. ", edited)) %>%
  mutate(edited = gsub("\\()", "", edited)) %>%
  mutate(edited = gsub("\\(\\(", "(", edited)) %>%
  mutate(edited = gsub("\\)\\)", ")", edited)) %>%
  mutate(edited = gsub(" x\\)", ")", edited)) %>%
  mutate(edited = gsub(" × [a-zA-Z]\\. × ", " × ", edited)) %>%
  mutate(edited = str_trim(str_squish(edited))) %>%
  filter(!(edited %in%
    c(
      "× Triticosecale Witt.", "× Triticosecale Wittmack", "Xanthosoma SD × Valentine",
      "× Triticosecale Wittm.", "× Ruttyruspolia × Phyllis van Heerden",
      "× Cardaminopsis HAYEK × Arabidopsis HEYNH.", "Bursera × kerbera hybrid with multijuga",
      "(Agropyron × Elymus ircutensis × Elymus) sp. Alt", "Vanda Josephine Van Brero × TMA",
      "Carex × Spec.", "Rosa cf. × Mill.", "Acer rubrum L. / Acer × freemanii E. Murray",
      "Zygopetalum × Rchb. f.", "Abelai × grandiflora", 
      "Q. faginea faginea × subsp. broteroi (Coutinho) A. Camus"
    ))) %>%
  mutate(edited = case_when(
    edited == "Aesculus × Carnea" ~ "Aesculus × carnea",
    edited == "AGATHOSMA CAPENSIS × ARIDA" ~ "Agathosma capensis × arida",
    edited == "Ambrosia dumosa × AChenopodifolia" ~ "Ambrosia dumosa × chenopodifolia",
    edited == "Amelanchier × Grandiflora" ~ "Amelanchier × grandiflora",
    edited == "Amelanchier × Wiegandii" ~ "Amelanchier × wiegandii",
    edited == "Annona × Atemoya" ~ "Annona × atemoya",
    edited == "Blechnum × Caudatum" ~ "Blechnum × caudatum",
    edited == "Blechnum × Confluens" ~ "Blechnum × confluens",
    edited == "Digitalis obscura × Lanata" ~ "Digitalis obscura × lanata",
    edited == "Dryopteris goldiana × Intermedia" ~ "Dryopteris goldiana × intermedia",
    edited == "Erigeron × Flahaultianum" ~ "Erigeron × flahaultianum",
    edited == "Galium verum × Mollugo" ~ "Galium verum × mollugo",
    edited == "Galium × Pomeranicum" ~ "Galium × pomeranicum",
    edited == "Impatiens × Pacifica" ~ "Impatiens × pacifica",
    edited == "Iris × Barthii" ~ "Iris × barthii",
    edited == "Mentha × VillosoNervata" ~ "Mentha × villosonervata",
    edited == "Musa × Paradisiaca" ~ "Musa × paradisiaca",
    edited == "Nicodemia madagascariensis × Buddleja × pikei" ~ 
      "Nicodemia madagascariensis × Buddleja pikei",
    edited == "Nicotiana × Sanderae" ~ "Nicotiana × sanderae",
    edited == "Potamogeton gramineus × Illinoensis" ~ "Potamogeton gramineus × illinoensis",
    edited == "Potamogeton × Fluitans" ~ "Potamogeton × fluitans",
    edited == "Potentilla arenaria × Tabernaemontani" ~ "Potentilla arenaria × tabernaemontani",
    edited == "Potentilla Crantzii × Tabernaemontani" ~ "Potentilla crantzii × tabernaemontani",
    edited == "Potentilla incana × Tabernaemontani" ~ "Potentilla incana × tabernaemontani",
    edited == "Quercus × Macdonaldii" ~ "Quercus × macdonaldii",
    edited == "Tilia × Vulgaris" ~ "Tilia × vulgaris",
    edited == "Inula × Yosezatoana Makino" ~ "Inula × yosezatoana Makino",
    edited == "Viola × dubia Wiesb. V. reichenbachiana × V. Riviniana" ~ "Viola × dubia Wiesb.",
    edited == "Viola × multicaulis Jord. (V. alba Bess. × odorata L.)" ~ 
      "Viola × multicaulis Jord.",
    edited == "× Cistus loreti Rouy et Foucaud C. ladanifero-monspeliensis" ~ 
      "× Cistus loreti Rouy et Foucaud × C. ladanifero-monspeliensis",
    edited == "× Senecio telonense Albert S. jacobeae × cineraria" ~ "Senecio telonense Albert",
    edited == 
      "(Miltoniopsis roezlii × M. vexillaria) × M. vexillaria) (M. vexillaria × M. roezlii)" ~ 
      "Miltoniopsis roezlii × M. vexillaria",
    edited == "(Elymus lanceolatus × Elymus glaucus) Elymus trachycaulus" ~ 
      "Elymus lanceolatus × Elymus glaucus",
    edited == "vriesea aff. × ruby" ~ "Vriesea aff. × ruby",
    TRUE ~ edited
  )) %>%
  mutate(edited = gsub("XTriticosecale|xTriticosecale", "× Triticosecale", edited)) %>%
  mutate(edited = gsub("^xMiscanthus", "Miscanthus", edited)) %>%
  mutate(parent = str_split(edited, pattern = "×")) %>%
  unnest("parent") %>%
  mutate(across(parent, ~ str_trim(str_squish(.)))) %>%
  mutate(parent = gsub("^\\(|\\)$", "", parent)) %>%
  mutate(across(parent, ~ str_trim(str_squish(.)))) %>%
  filter(parent != "") %>%
  mutate(parent = gsub(" [a-zA-Z]$", "", parent)) %>%
  mutate(across(parent, ~ str_trim(str_squish(.)))) %>%
  rowwise() %>%
  mutate(genus = str_split(
    gsub("^× ", "", edited),
    pattern = " ", simplify = TRUE
  )[1]) %>%
  mutate(genus = str_to_sentence(genus)) %>%
  filter(parent != genus) %>%
  mutate(genus = gsub("^\\(", "", genus)) %>%
  mutate(parent = gsub("^[A-Za-z]\\.|^[A-Za-z] ", paste(genus, " "), parent)) %>%
  mutate(across(parent, ~ str_trim(str_squish(.)))) %>%
  mutate(parent = ifelse(
    grepl("^[a-z]", parent) & !grepl("^subsp\\.", parent),
    paste(genus, parent), parent
  )) %>%
  mutate(parent = ifelse(grepl("subsp\\. ×", edited), NA, parent)) %>%
  filter(str_count(parent, pattern = " ") > 0) %>%
  mutate(parent = ifelse(gsub(" × ", " ", edited) == parent, NA, parent)) %>%
  mutate(parent = ifelse(grepl(" aff\\. ×", " ", edited), NA, parent)) %>%
  mutate(parent = ifelse(
    grepl("^subsp\\.", parent),
    paste(
      str_split(gsub("^× ", "", edited), pattern = " ", simplify = TRUE)[1],
      str_split(gsub("^× ", "", edited), pattern = " ", simplify = TRUE)[2], parent
    ), parent
  )) %>%
  select(-genus) %>%
  ungroup() %>%
  mutate(across(parent, ~ str_trim(str_squish(.)))) %>%
  distinct()


# combine -----------------------------------------------------------------------------------------------

hybrids_gbif_occ_all <- bind_rows(
  hybrids_gbif_occ1, hybrids_gbif_occ2, hybrids_gbif_occ3
) %>%
  distinct()


# match to GBIF -----------------------------------------------------------------------------------------

source(here("r_functions/fn_gbif_matchTaxa.R"))
gbif.cmpfn <- compiler::cmpfun(gbif.fn)
gbif_matched <- gbif.cmpfn(
  taxonName = unique(hybrids_gbif_occ_all$edited),
  taxonID = seq_along(unique(hybrids_gbif_occ_all$edited))
) %>%
  select(taxonName, confidence, matchType, scientificName, usageKey) %>%
  rename(edited = taxonName) %>%
  distinct()

# names with multiple IDs
gbif_matched %>%
  select(scientificName, usageKey) %>%
  distinct() %>%
  count(.$scientificName) %>%
  filter(n > 1)

# exclude potentially incorrect matches
gbif_matched %<>%
  filter(!(
    grepl(" subsp\\.", scientificName) &
      !grepl("subsp\\.|var\\.", edited) &
      matchType == "FUZZY"))
gbif_matched %<>%
  filter(!(
    grepl(" var\\.", scientificName) &
      !grepl("×", scientificName) &
      !grepl("var\\.|subsp\\.", edited) &
      matchType == "FUZZY"))

# match parent taxa
parent_gbif_matched <- gbif.cmpfn(
  taxonName = unique(hybrids_gbif_occ_all$parent),
  taxonID = seq_along(unique(hybrids_gbif_occ_all$parent))
) %>%
  select(taxonName, confidence, matchType, scientificName, usageKey) %>%
  rename_with(~ paste0("parent_", .), .cols = everything()) %>%
  rename(parent = parent_taxonName) %>%
  distinct()

# names with multiple IDs
parent_gbif_matched %>%
  select(parent_scientificName, parent_usageKey) %>%
  distinct() %>%
  count(.$parent_scientificName) %>%
  filter(n > 1)

hybrids_gbif_occ_all_matched <- hybrids_gbif_occ_all %>%
  left_join(gbif_matched) %>%
  left_join(parent_gbif_matched) %>%
  filter(!is.na(usageKey)) %>%
  select(-edited, -parent) %>%
  distinct()

hybrids_gbif_occ_all %<>%
  left_join(gbif_matched) %>%
  left_join(parent_gbif_matched) %>%
  filter(is.na(usageKey)) %>%
  bind_rows(hybrids_gbif_occ_all_matched) %>%
  write_csv2(here("data/hybrid_names_from_gbif_occurrences.csv"))

# hybrids_w_parents <- hybrids_gbif_occ_all %>%
#   filter(!is.na(parent_scientificName))
# hybrids_wo_parents <- hybrids_gbif_occ_all %>%
#   filter(is.na(parent_scientificName)) %>%
#   filter(!(scientificName %in% hybrids_w_parents$scientificName))
# hybrids_gbif_occ_all <- bind_rows(hybrids_w_parents, hybrids_wo_parents)
