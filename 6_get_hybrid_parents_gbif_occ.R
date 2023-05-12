# Marina Golivets
# BiCIKL Hackathon (Sept 20-24, 2021, Meise, Belgium)
# Clean hybrid names originating from GBIF occurrence data
# Extract parent taxa


# The list was split into three parts:
# 1) list of hybrid names + hybrid formulas (with "=")
# 2) list of hybrid names
# 3) the rest


# part 1 -------------------------------------------------------------------------------------------
hybrids_gbif_occ1 <- tibble(
  verbatim = read_lines(
    here("input_datasets/hybridnamesWithHybridFormulasFromOccurences.txt")
  )
) %>%
  mutate(
    edited = str_remove_all(verbatim, "\\?|- Festul asc") %>%
      str_replace_all(" ssp\\. ", " subsp. ") %>%
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
  mutate(across(everything(), ~ str_trim(str_squish(.))))

# part 2 -------------------------------------------------------------------------------------------
hybrids_gbif_occ2 <- tibble(
  verbatim = read_lines(
    here("input_datasets/hybridnamesFromOccurences.txt")
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
  )) %>%
  filter(
    str_detect(edited, " × [A-Z][a-z]", negate = TRUE),
    str_detect(edited, "× sp\\.$|× hybr\\.$|× L\\.$", negate = TRUE)
  )

# part 3 -------------------------------------------------------------------------------------------
hybrids_gbif_occ3 <- data.table::fread(
  here("input_datasets/2021-09-23_14-47_hybridnamesUniq2_parents.csv")
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
        ".*or possibly", "="),
        collapse = "|"
      ),
      negate = TRUE
    ),
  ) %>%
  mutate(
    edited = edited %>%
      str_replace_all("× ×|× (?i)x | (?i)x (?i)x ", "×") %>%
      str_replace_all(" (?i)x |^(?i)x |×|x_|^(?i)x-", " × ") %>%
      str_remove_all("[0-9]|#|;|&#|N\\/A") %>%
      str_squish() %>%
      str_trim() %>%
      str_remove_all(
        str_c(
         c( "^[a-zA-Z]+ceae ", "^[A-Z]+CEAE ",
          "Complex hybrid", "doubled", "sterile triploid", "possible hybrid",
          "autotetraploid", "diploid", "\\(allo", "suggested hybrid", "form A", "form B",
          "unranked", ".*or possibly"),
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
          "<ALL>", "× [a-zA-Z]$", "× [a-zA-Z]\\.$", "× ssp\\.$"),
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
      " nssp\\." = " nothosubsp. ",
      " subvar\\. × " = " subvar\\. ",
      "Zostera nana subvar. X Zostera marina" = "Zostera nana × Zostera marina",
      "-\\)" = ")"
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
  mutate(
    edited = edited %>%
      str_replace_all(
        c(
          "XTriticosecale|xTriticosecale" = "× Triticosecale",
          "^xMiscanthus" = "Miscanthus"
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
hybrids_gbif_occ_all <- bind_rows(
  hybrids_gbif_occ1, hybrids_gbif_occ2, hybrids_gbif_occ3
) %>%
  distinct() %>%
  # group_by(verbatim) %>%
  # filter(n() > 1) %>%
  write_csv2(here("output_datasets/hybrids_gbif_occurrences.csv"))