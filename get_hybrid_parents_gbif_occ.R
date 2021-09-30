# Marina Golivets
# BiCIKL Hackathon (Sept 20-24, 2021, Meise)
# Clean hybrid names originating from GBIF occurrence data
# Extract parent taxa


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
  mutate(edited = gsub("×| x |^x ", " × ", edited, ignore.case = TRUE)) %>%
  mutate(across(edited, ~ str_trim(str_squish(.)))) %>%
  mutate(edited = gsub("Erigeron xhuelsenii", "Erigeron × huelsenii", edited)) %>%
  rowwise() %>%
  mutate(genus = str_split(edited, pattern = " ", simplify = TRUE)[1]) %>%
  separate(edited, into = c("a", "b"), sep = "=|<=>", remove = FALSE) %>%
  mutate(across(a:b, ~ str_trim(str_squish(.)))) %>%
  mutate(a = gsub("\\($|\\[$", "", a)) %>%
  mutate(b = gsub("^\\(|\\)$|\\]$", "", b)) %>%
  mutate(across(a:b, ~ str_trim(str_squish(.)))) %>%
  mutate(b = gsub("c. dioeca × davalliana) Ascherson & Graebner", 
                  "C. dioeca × davalliana", b)) %>%
  mutate(b = ifelse(
    grepl("^[a-zA-Z]\\.|^× [a-zA-Z]\\.|^× [a-z]+ |^[A-Z] |^[a-zA-Z]{2}\\.", b),
    paste(genus, gsub("^[a-zA-Z]\\.|^× [a-zA-Z]\\.|^× [a-z]+ |^[A-Z] |^[a-zA-Z]{2}\\.", "", b)), b
  )) %>%
  rowwise() %>%
  mutate(b = ifelse(grepl("^[a-z]", b), paste(genus, b), b)) %>%
  mutate(b = gsub("^× ", "", b)) %>%
  mutate(across(a:b, ~ str_trim(str_squish(.)))) %>%
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
  mutate(formula = ifelse(is.na(formula) & grepl("×", a), a, formula)) %>%
  mutate(formula = case_when(
    !(formula %in% c("Platanus × hybrida", "Rumex × pratensis M. et K.")) ~ formula)
    ) %>%
  select(verbatim, a, b, formula) %>%
  pivot_longer(a:b, values_to = "edited") %>%
  select(-name) %>%
  filter(!is.na(edited)) %>%
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

hybrids_gbif_occ3 <- data.table::fread(here("data/2021-09-23_14-47_hybridnamesUniq2_parents.csv")) %>%
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
  mutate(edited = gsub(" subvar\\. × ", " subvar\\. ", edited)) %>%
  mutate(edited = case_when(
    verbatim == "Zostera nana subvar. X Zostera marina" ~ "Zostera nana × Zostera marina",
    TRUE ~ edited
    )) %>%
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


# combine ------------------------------------------------------------------------------------------

hybrids_gbif_occ_all <- bind_rows(
  hybrids_gbif_occ1, hybrids_gbif_occ2, hybrids_gbif_occ3
) %>%
  distinct() %>%
  write_csv2(here("data/hybrids_gbif_occurrences.csv"))