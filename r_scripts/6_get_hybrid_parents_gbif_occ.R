# Marina Golivets
# BiCIKL Hackathon (Sept 20-24, 2021, Meise, Belgium)
# Clean hybrid names originating from GBIF occurrence data
# Extract parent taxa

# The list was split into three parts:
# 1) list of hybrid names + hybrid formulas (with "=")
# 2) list of hybrid names
# 3) the rest

# part 1 -------------------------------------------------------------------------------------------
hybr_gbif_occ1 <- tibble(
  verbatim = read_lines(
    here("input/hybridnamesWithHybridFormulasFromOccurences.txt")
  )
) %>%
  mutate(
    edited = str_remove_all(verbatim, "\\?|- Festul asc") %>%
      str_replace_all(" ssp\\.", " subsp. ") %>%
      str_replace_all("\\.x |\\.x\\.", ". × ") %>%
      str_replace_all("×| (?i)x |^(?i)x ", " × ") %>%
      str_squish() %>%
      str_trim() %>%
      str_replace_all("Erigeron xhuelsenii", "Erigeron × huelsenii")
  ) %>%
  mutate(genus = str_extract(edited, "\\w*")) %>%
  separate(edited, into = c("a", "b"), sep = "=|<=>", remove = FALSE) %>%
  mutate(
    across(a:b, ~ str_trim(str_squish(.))),
    a = str_remove_all(a, "\\($|\\[$"),
    b = str_remove_all(b, "^\\(|\\)$|\\]$"),
    across(a:b, ~ str_trim(str_squish(.)))
  ) %>%
  mutate(
    b = str_replace(
      b, "c. dioeca × davalliana\\) Ascherson & Graebner",
      "C. dioeca × davalliana"
    ),
    b = if_else(
      str_detect(b, "^[a-zA-Z]\\.|^× [a-zA-Z]\\.|^× [a-z]+ |^[A-Z] |^[a-zA-Z]{2}\\."),
      paste(
        genus,
        str_remove_all(b, "^[a-zA-Z]\\.|^× [a-zA-Z]\\.|^× [a-z]+ |^[A-Z] |^[a-zA-Z]{2}\\.")
      ),
      b
    ),
    b = str_detect(b, "^[a-z]") %>%
      if_else(paste(genus, b), b) %>%
      str_remove("^× ")
  ) %>%
  mutate(
    across(a:b, ~ str_trim(str_squish(.))),
    a_f = str_detect(a, paste(genus, "×")),
    b_f = str_detect(b, paste(genus, "×")),
    c_f = str_detect(a, "×"),
    d_f = str_detect(b, "×")
  ) %>%
  mutate(formula = case_when(
    a_f == FALSE & b_f == TRUE ~ a,
    b_f == FALSE & a_f == TRUE ~ b,
    b_f == FALSE & a_f == FALSE & c_f == TRUE ~ a,
    b_f == FALSE & a_f == FALSE & d_f == TRUE ~ b
  )) %>%
  mutate(formula = ifelse(is.na(formula) & grepl("×", a), a, formula)) %>%
  mutate(formula = case_when(
    !(formula %in% c("Platanus × hybrida", "Rumex × pratensis M. et K.")) ~ formula
  )) %>%
  select(verbatim, a, b, formula) %>%
  pivot_longer(a:b, values_to = "edited") %>%
  select(-name) %>%
  filter(!is.na(edited)) %>%
  mutate(parent = str_split(formula, pattern = "×")) %>%
  unnest("parent") %>%
  mutate(across(parent, ~ str_trim(str_squish(.)))) %>%
  mutate(
    genus = str_extract(formula, "\\w*"),
    parent = case_when(
      grepl("^[a-z]", parent) ~ paste(genus, parent),
      grepl("^[A-Z]\\.", parent) ~
        paste(genus, str_remove(parent, "^[A-Z]\\.")),
      TRUE ~ parent
    )
  ) %>%
  select(verbatim, edited, parent) %>%
  mutate(across(everything(), ~ str_squish(.) %>% str_trim()))

# part 2 -------------------------------------------------------------------------------------------
hybr_gbif_occ2 <- tibble(
  verbatim = read_lines(
    here("input/hybridnamesFromOccurences.txt")
  )
) %>%
  filter(
    str_detect(verbatim, "'", negate = TRUE),
    !(verbatim %in% c(
      "S.ilicifolius ×S.soldanelloides",
      "S.multisetus ×S.triangulatus",
      "pelargonium ×l.h.bailey",
      "L.pilosa ×L.cuneata",
      "Rumex ×longif.-pseudonatr.",
      "Abelai ×grandiflora"
    ))
  ) %>%
  mutate(
    edited = verbatim %>%
      str_replace_all("×", " × ") %>%
      str_squish() %>%
      str_trim()
  ) %>%
  mutate(edited = case_when(
    verbatim == "Tripsacum ×maizar-zopilotense-hybrid)" ~ "Tripsacum maizar × zopilotense",
    TRUE ~ edited
  ) %>%
    str_replace("Salicaceae", "Salix")
  ) %>%
  filter(
    str_detect(edited, " × [A-Z][a-z]", negate = TRUE),
    str_detect(edited, "× sp\\.$|× hybr\\.$|× L\\.$", negate = TRUE)
  )

# part 3 -------------------------------------------------------------------------------------------
hybr_gbif_occ3 <- data.table::fread(
  here("input/2021-09-23_14-47_hybridnamesUniq2_parents.csv")
) %>%
  transmute(
    verbatim = input_strings,
    edited = str_squish(verbatim) %>%
      str_trim()
  ) %>%
  filter(
    edited %>% str_detect(
      str_c(
        c(" '", "' ", "´", "‘", "’", " cv\\.", " f\\. [A-Z]", " var\\. [A-Z] ", " or ",
        ".*or possibly", "=", "\\) .* \\(", ">"),
        collapse = "|"
      ),
      negate = TRUE
    ),
  ) %>%
  mutate(
    edited = edited %>%
      str_replace_all(" x x |× ×|× (?i)x | (?i)x (?i)x | x\\. | × [A-Z]\\. ", "×") %>%
      str_replace_all(" (?i)x |^(?i)x |×|x_|^(?i)x-", " × ") %>%
      str_remove_all("[0-9]|#|;|&#|N\\/A") %>%
      str_squish() %>%
      str_trim() %>%
      str_remove_all(
        str_c(
         c( "^[a-zA-Z]+ceae ", "^[A-Z]+CEAE ",
          "Complex hybrid", "doubled", "sterile triploid", "possible hybrid",
          "autotetraploid", "diploid", "\\(allo", "suggested hybrid", "form A", "form B",
          "unranked", ".*or possibly", " x$"),
          collapse = "|"
        )
      ) %>%
      str_replace_all("Cliffortia sp\\.", "Cliffortia")
  ) %>%
  filter(
    edited %>%
      str_detect(
        str_c(
          c("× \\?$", " × [A-Z][a-z]+ sp\\.", "× unknown$", " unknown species", "indet\\. sp\\.",
          "<ALL>", "× [a-zA-Z]$", "× [a-zA-Z]\\.$", "× ssp\\.$", "× var\\.",
          "^[A-Z]. ", "Polystichum × Sa.Kuratae", "Loch Valley"),
          collapse = "|"
        ),
        negate = TRUE
      )
  ) %>% # exclude hybrid formulas with partially unknown parents
  mutate(
    edited = str_remove(edited, "unknown") %>%
      str_remove_all("\\?|_|\\[|\\]| -|- |\\+|-$|--$|pro sp\\.|pro\\. sp\\.")
  ) %>%
  filter(str_detect(edited, " sp\\.", negate = TRUE)) %>%
  mutate(edited = str_replace_all(
    edited,
    c(
      " ssp\\." = " subsp. ",
      " ssp " = " subsp. ",
      " nssp\\." = " nothosubsp. ",
      " subvar\\. × " = " subvar\\. ",
      "Zostera nana subvar. X Zostera marina" = "Zostera nana × Zostera marina",
      "-\\)" = ")",
      "^abelia" = "Abelia",
      "Prunus \\(Prunus" = "Prunus",
      " abies " = " Abies "
    )
  ) %>%
    str_squish() %>%
    str_trim() %>%
    str_remove_all("\\(\\)") %>%
    str_replace_all(
      c(
        "\\(\\(" = "(",
        "\\)\\)" = ")",
        "\\)\\)" = ")",
        " x\\)" = ")",
        " × [a-zA-Z]\\. × " = " × "
      )
    ) %>%
    str_squish() %>%
    str_trim()) %>%
  filter(!(edited %in%
    c(
      "× Triticosecale Witt.", "× Triticosecale Wittmack", "Xanthosoma SD × Valentine",
      "× Triticosecale Wittm.", "× Ruttyruspolia × Phyllis van Heerden",
      "× Cardaminopsis HAYEK × Arabidopsis HEYNH.", "Bursera × kerbera hybrid with multijuga",
      "(Agropyron × Elymus ircutensis × Elymus) sp. Alt", "Vanda Josephine Van Brero × TMA",
      "Carex × Spec.", "Rosa cf. × Mill.", "Acer rubrum L. / Acer × freemanii E. Murray",
      "Zygopetalum × Rchb. f.", "Abelai × grandiflora",
      "Q. faginea faginea × subsp. broteroi (Coutinho) A. Camus",
      "pelargonium × l.h.bailey", "Genus × spec.", "Hemerocallis × Stella de Oro",
      "Citrus × (L.) Osbeck", "Angustifolia  × L. Latifolia", "Vriesea × E. Morren",
      "Elymus × [A-Z]", "Rosa canina L. s. l. ( × )", "Vernonia platensis var. × (Vell.) Rusby"
    ))) %>%
  mutate(edited = case_match(edited,
    "Aesculus × Carnea" ~ "Aesculus × carnea",
    "AGATHOSMA CAPENSIS × ARIDA" ~ "Agathosma capensis × arida",
    "Ambrosia dumosa × AChenopodifolia" ~ "Ambrosia dumosa × chenopodifolia",
    "Amelanchier × Grandiflora" ~ "Amelanchier × grandiflora",
    "Amelanchier × Wiegandii" ~ "Amelanchier × wiegandii",
    "Annona × Atemoya" ~ "Annona × atemoya",
    "Blechnum × Caudatum" ~ "Blechnum × caudatum",
    "Blechnum × Confluens" ~ "Blechnum × confluens",
    "Digitalis obscura × Lanata" ~ "Digitalis obscura × lanata",
    "Dryopteris goldiana × Intermedia" ~ "Dryopteris goldiana × intermedia",
    "Erigeron × Flahaultianum" ~ "Erigeron × flahaultianum",
    "Galium verum × Mollugo" ~ "Galium verum × mollugo",
    "Galium × Pomeranicum" ~ "Galium × pomeranicum",
    "Impatiens × Pacifica" ~ "Impatiens × pacifica",
    "Iris × Barthii" ~ "Iris × barthii",
    "Mentha × VillosoNervata" ~ "Mentha × villosonervata",
    "Musa × Paradisiaca" ~ "Musa × paradisiaca",
    "Nicodemia madagascariensis × Buddleja × pikei" ~
      "Nicodemia madagascariensis × Buddleja pikei",
    "Nicotiana × Sanderae" ~ "Nicotiana × sanderae",
    "Potamogeton gramineus × Illinoensis" ~ "Potamogeton gramineus × illinoensis",
    "Potamogeton × Fluitans" ~ "Potamogeton × fluitans",
    "Potentilla arenaria × Tabernaemontani" ~ "Potentilla arenaria × tabernaemontani",
    "Potentilla Crantzii × Tabernaemontani" ~ "Potentilla crantzii × tabernaemontani",
    "Potentilla incana × Tabernaemontani" ~ "Potentilla incana × tabernaemontani",
    "Quercus × Macdonaldii" ~ "Quercus × macdonaldii",
    "Tilia × Vulgaris" ~ "Tilia × vulgaris",
    "Inula × Yosezatoana Makino" ~ "Inula × yosezatoana Makino",
    "Viola × dubia Wiesb. V. reichenbachiana × V. Riviniana" ~ "Viola × dubia Wiesb.",
    "Viola × multicaulis Jord. (V. alba Bess. × odorata L.)" ~
      "Viola × multicaulis Jord.",
    "× Cistus loreti Rouy et Foucaud C. ladanifero-monspeliensis" ~
      "× Cistus loreti Rouy et Foucaud × C. ladanifero-monspeliensis",
    "× Senecio telonense Albert S. jacobeae × cineraria" ~ "Senecio telonense Albert",
      "(Miltoniopsis roezlii × M. vexillaria) × M. vexillaria) (M. vexillaria × M. roezlii)" ~
      "Miltoniopsis roezlii × M. vexillaria",
    "(Elymus lanceolatus × Elymus glaucus) Elymus trachycaulus" ~
      "Elymus lanceolatus × Elymus glaucus",
    "vriesea aff. × ruby" ~ "Vriesea aff. × ruby",
    "Elymus elymus × littoreus" ~ "Elymus × littoreus",
    "Elymus × littoreus (C.F. Schumach.) Lambinon in D" ~ "Elymus × littoreus (C.F. Schumach.) Lambinon",
    "Lobelia × L. cardinalis L." ~ "Lobelia × cardinalis L.",
    "Quercus × Q. ithaburensis subsp. macrolepis" ~ "Quercus × ithaburensis subsp. macrolepis",
    "Amaranthus × A. tucsonensis Henrickson" ~ "Amaranthus × tucsonensis Henrickson",
    "Cardamine × C. maxima" ~ "Cardamine × maxima",
    "Polystichum × Dbicknellii (Christ) Hahne" ~ "Polystichum × bicknellii (Christ) Hahne",
    "Platanthera × Canbyi (Ames) Luer"  ~"Platanthera × canbyi (Ames) Luer"  ,
    "Pistacia ×Saportae Buznat" ~"Pistacia ×saportae Buznat" ,
    "Nuphar × Nuphar rubrodiscum" ~ "Nuphar × rubrodiscum",
    "Medicago Martyn × varia" ~ "Medicago × varia",
    "vriesea × brueggemannii matos & crespo" ~ "Vriesea ×brueggemannii Matos & Crespo",
    "tilia × moltkei L." ~ "Tilia × moltkei L.",
    "prunus × yedoensis matsum." ~ "Prunus ×yedoensis Matsum.",
    "Acacia hugelii x a. glauca" ~ "Acacia hugelii × A. glauca",
    .default = edited
  )) %>%
  mutate(
    edited = edited %>%
      str_replace_all(
        c(
          "XTriticosecale|xTriticosecale" = "× Triticosecale",
          "^xMiscanthus" = "Miscanthus",
          "Pelargonium × P." = "Pelargonium ×"
        )
      ),
    parent = str_split(edited, pattern = "×")
  ) %>%
  unnest("parent") %>%
  mutate(across(parent, ~ str_trim(str_squish(.)))) %>%
  mutate(
    parent = parent %>%
      str_remove_all("^\\(|\\)$") %>%
      str_squish() %>%
      str_trim()
  ) %>%
  filter(parent != "") %>%
  mutate(
    parent = parent %>%
      str_remove_all(" [a-zA-Z]$") %>%
      str_squish() %>%
      str_trim(),
    genus = edited %>%
      str_remove("^× ") %>%
      str_extract("\\w*") %>%
      str_to_sentence()
  ) %>%
  filter(parent != genus) %>%
  mutate(
    genus = str_remove(genus, "^\\("),
    parent = str_replace(parent, "^[A-Za-z]\\.|^[A-Za-z] ", paste(genus, " ")) %>%
      str_squish() %>%
      str_trim() %>%
      if_else(
        str_detect(., "^[a-z]") & str_detect(., "^subsp\\.", negate = TRUE),
        paste(genus, .), .
      )
  ) %>%
  mutate(parent = if_else(str_detect(edited, "subsp\\. ×"), NA_character_, parent)) %>%
  filter(str_count(parent, pattern = " ") > 0) %>%
  mutate(
    parent = case_when(
      str_replace_all(edited, " × ", " ") == parent ~ NA_character_,
      str_replace_all(edited, " aff\\. ×", " ") == parent ~ NA_character_,
      TRUE ~ parent
    ) %>%
      if_else(
        str_detect(., "^subsp\\."),
        paste(
          edited %>%
            str_remove("^× ") %>%
            str_extract("(\\w+\\s+\\w+)"),
          .
        ), .
      ) %>%
      str_squish() %>%
      str_trim()
  ) %>%
  select(-genus) %>%
  distinct()

# combine ------------------------------------------------------------------------------------------
hybr_gbif_occ_all1 <- bind_rows(
  hybr_gbif_occ1, hybr_gbif_occ2, hybr_gbif_occ3
) %>%
  mutate(
    type = if_else(str_detect(edited, "^[A-Za-z]+ × |^× "), "name", "formula") %>%
      if_else(str_detect(edited, "^[A-Za-z]+ × .*×"), "formula", .) %>%
      if_else(str_detect(edited, "×", negate = TRUE), "name", .),
    edited = if_else(
      type == "name", str_replace_all(edited, c(" × " = " ×", "^× " = "×")), edited) %>%
      str_replace(" × × ", " × "),
    source = "gbif occurrences",
    dist = stringdist::stringdist(edited, parent, method = "lv"),
    parent = ifelse(dist > 3, parent, NA_character_)
    ) %>%
  select(-dist) %>%
  distinct()

names <- hybr_gbif_occ_all1 %>% 
  filter(type == "name")
formulas <- hybr_gbif_occ_all1 %>% 
  filter(type == "formula", !(verbatim %in% names$verbatim))
hybr_gbif_occ_all2 <- bind_rows(names, formulas) %>%
  mutate(edited = if_else(type == "name", str_replace(edited, " × ", " ×"), edited)) %>%
  write_csv2(here("output/hybrids_gbif_occurrences.csv"))
