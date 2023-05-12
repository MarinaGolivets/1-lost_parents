# Marina Golivets
# BiCIKL Hackathon
# Extract hybrid taxon names from various sources


devtools::install_github("ropensci/chromer")
# parse the list from Quentin ---------------------------------------------------------------------------


all_hybrids2 <- all_hybrids1 %>%
  bind_rows(hybrids_gbif_occ1) %>%
  mutate(verbatimTaxonName = gsub("× [A-Z]\\.", "× ", verbatimTaxonName)) %>%
  mutate(parent_verbatimTaxonName = ifelse(
    parent_verbatimTaxonName == verbatimTaxonName, NA, parent_verbatimTaxonName
  )) %>%
  mutate(parent_verbatimTaxonName = ifelse(
    gsub(" × ", " ", verbatimTaxonName) == parent_verbatimTaxonName, NA, parent_verbatimTaxonName
  )) %>%
  distinct()




all_hybrids3 <- all_hybrids2 %>%
  bind_rows(hybrids_ecoflora) %>%
  mutate(verbatimTaxonName = gsub("× [A-Z]\\.", "× ", verbatimTaxonName)) %>%
  mutate(parent_verbatimTaxonName = ifelse(
    parent_verbatimTaxonName == verbatimTaxonName, NA, parent_verbatimTaxonName
  )) %>%
  mutate(parent_verbatimTaxonName = ifelse(
    gsub(" × ", " ", verbatimTaxonName) == parent_verbatimTaxonName, NA, parent_verbatimTaxonName
  )) %>%
  distinct() %>%


all_hybrids4 <- all_hybrids3 %>%
  bind_rows(hybrids_ipni %>% select(verbatimTaxonName, parent_verbatimTaxonName)) %>%
  mutate(parent_verbatimTaxonName = ifelse(parent_verbatimTaxonName == "", NA, parent_verbatimTaxonName)) %>%
  distinct()


# parse hybrid names from GBIF occurrences --------------------------------------------------------------
all_hybrids5 <- all_hybrids4 %>%
  bind_rows(hybrids_gbif_occ2) %>%
  mutate(verbatimTaxonName = gsub("× [A-Z]\\.", "× ", verbatimTaxonName)) %>%
  mutate(parent_verbatimTaxonName = ifelse(
    parent_verbatimTaxonName == verbatimTaxonName, NA, parent_verbatimTaxonName
  )) %>%
  mutate(parent_verbatimTaxonName = ifelse(
    gsub(" × ", " ", verbatimTaxonName) == parent_verbatimTaxonName, NA, parent_verbatimTaxonName
  )) %>%
  distinct()




all_hybrids6 <- all_hybrids5 %>%
  bind_rows(gbif_chl) %>%
  mutate(parent_verbatimTaxonName = ifelse(
    parent_verbatimTaxonName == verbatimTaxonName, NA, parent_verbatimTaxonName
  )) %>%
  mutate(parent_verbatimTaxonName = ifelse(
    gsub(" × ", " ", verbatimTaxonName) == parent_verbatimTaxonName, NA, parent_verbatimTaxonName
  )) %>%
  mutate(parent_verbatimTaxonName = gsub("ssp\\.", "subsp.", parent_verbatimTaxonName)) %>%
  mutate(verbatimTaxonName = gsub("ssp\\.", "subsp.", verbatimTaxonName)) %>%
  distinct()

hybrids_w_parents <- all_hybrids6 %>%
  filter(!is.na(parent_verbatimTaxonName))
hybrids_wo_parents <- all_hybrids6 %>%
  filter(is.na(parent_verbatimTaxonName)) %>%
  filter(!(verbatimTaxonName %in% hybrids_w_parents$verbatimTaxonName))
all_hybrids6 <- bind_rows(hybrids_w_parents, hybrids_wo_parents)


# match hybrid taxon names to GBIF backbone --------------------------------------------
source(here("r_functions/fn_gbif_matchTaxa.R"))
gbif.cmpfn <- compiler::cmpfun(gbif.fn)
gbif_matched <- gbif.cmpfn(
  taxonName = unique(all_hybrids6$verbatimTaxonName),
  taxonID = seq_along(unique(all_hybrids6$verbatimTaxonName))
) %>%
  select(taxonName, confidence, matchType, scientificName, usageKey) %>%
  rename(verbatimTaxonName = taxonName) %>%
  distinct()

# gbif_matched %>%
#   select(scientificName, usageKey) %>%
#   distinct() %>%
#   count(.$scientificName) %>%
#   filter(n > 1)

# exclude incorrect matches
gbif_matched %<>%
  filter(!(
    grepl(" subsp\\.", scientificName) &
      !grepl("subsp\\.|var\\.", verbatimTaxonName) &
      matchType == "FUZZY"))

gbif_matched %<>%
  filter(!(
    grepl(" var\\.", scientificName) &
      !grepl("×", scientificName) &
      !grepl("var\\.|subsp\\.", verbatimTaxonName) &
      matchType == "FUZZY"))

parent_gbif_matched <- gbif.cmpfn(
  taxonName = unique(all_hybrids6$parent_verbatimTaxonName),
  taxonID = seq_along(unique(all_hybrids6$parent_verbatimTaxonName))
) %>%
  select(taxonName, confidence, matchType, scientificName, usageKey) %>%
  rename(verbatimTaxonName = taxonName) %>%
  rename_with(~ paste0("parent_", .), .cols = everything()) %>%
  distinct()

# parent_gbif_matched %>%
#   select(parent_scientificName, parent_usageKey) %>%
#   distinct() %>%
#   count(.$parent_scientificName) %>%
#   filter(n > 1)

all_hybrids6_matched <- all_hybrids6 %>%
  left_join(gbif_matched) %>%
  left_join(parent_gbif_matched) %>%
  filter(!is.na(usageKey)) %>%
  select(-verbatimTaxonName, -parent_verbatimTaxonName) %>%
  distinct()

hybrids_w_parents <- all_hybrids6_matched %>%
  filter(!is.na(parent_scientificName))
hybrids_wo_parents <- all_hybrids6_matched %>%
  filter(is.na(parent_scientificName)) %>%
  filter(!(scientificName %in% hybrids_w_parents$scientificName))
all_hybrids6_matched <- bind_rows(hybrids_w_parents, hybrids_wo_parents)

all_hybrids6 %<>%
  left_join(gbif_matched) %>%
  left_join(parent_gbif_matched) %>%
  filter(is.na(usageKey)) %>%
  mutate(scientificName = verbatimTaxonName) %>%
  bind_rows(all_hybrids6_matched)

hybrids_w_parents <- all_hybrids6 %>%
  filter(!is.na(parent_scientificName))
hybrids_wo_parents <- all_hybrids6 %>%
  filter(is.na(parent_scientificName)) %>%
  filter(!(scientificName %in% hybrids_w_parents$scientificName))
all_hybrids6 <- bind_rows(hybrids_w_parents, hybrids_wo_parents)

all_hybrids6 %<>%
  mutate(parent_scientificName = ifelse(
    is.na(parent_scientificName), parent_verbatimTaxonName, parent_scientificName
  )) %>%
  distinct()

# get hybrid names from Wikidata ---------------------------------------------------------------


all_hybrids7 <- all_hybrids6 %>%
  bind_rows(hybrids_wiki) %>%
  select(-verbatimTaxonName, -parent_verbatimTaxonName) %>%
  distinct()


# get the rest of the names from GBIF occurrence data ---------------------------------------


gbif_occ3_matched <- gbif.cmpfn(
  taxonName = unique(hybrids_gbif_occ3$verbatimTaxonName),
  taxonID = seq_along(unique(hybrids_gbif_occ3$verbatimTaxonName))
) %>%
  select(taxonName, confidence, matchType, scientificName, usageKey) %>%
  rename(verbatimTaxonName = taxonName) %>%
  distinct()

parent_gbif_occ3_matched <- gbif.cmpfn(
  taxonName = unique(hybrids_gbif_occ3$parent_verbatimTaxonName),
  taxonID = seq_along(unique(hybrids_gbif_occ3$parent_verbatimTaxonName))
) %>%
  select(taxonName, confidence, matchType, scientificName, usageKey) %>%
  rename(verbatimTaxonName = taxonName) %>%
  rename_with(~ paste0("parent_", .), .cols = everything()) %>%
  distinct()

# exclude incorrect matches
gbif_occ3_matched %<>%
  filter(!(
    grepl(" subsp\\.", scientificName) &
      !grepl("subsp\\.|var\\.", verbatimTaxonName) &
      matchType == "FUZZY"))

gbif_occ3_matched %<>%
  filter(!(
    grepl(" var\\.", scientificName) &
      !grepl("×", scientificName) &
      !grepl("var\\.|subsp\\.", verbatimTaxonName) &
      matchType == "FUZZY"))

hybrids_gbif_occ3_matched <- hybrids_gbif_occ3 %>%
  left_join(gbif_occ3_matched) %>%
  left_join(parent_gbif_occ3_matched) %>%
  filter(!is.na(usageKey)) %>%
  select(-verbatimTaxonName, -parent_verbatimTaxonName) %>%
  distinct()

hybrids_w_parents <- hybrids_gbif_occ3_matched %>%
  filter(!is.na(parent_scientificName))
hybrids_wo_parents <- hybrids_gbif_occ3_matched %>%
  filter(is.na(parent_scientificName)) %>%
  filter(!(scientificName %in% hybrids_w_parents$scientificName))
hybrids_gbif_occ3_matched <- bind_rows(hybrids_w_parents, hybrids_wo_parents)

hybrids_gbif_occ3 %<>%
  left_join(gbif_occ3_matched) %>%
  left_join(parent_gbif_occ3_matched) %>%
  filter(is.na(usageKey)) %>%
  mutate(scientificName = verbatimTaxonName) %>%
  bind_rows(hybrids_gbif_occ3_matched)

hybrids_w_parents <- hybrids_gbif_occ3 %>%
  filter(!is.na(parent_scientificName))
hybrids_wo_parents <- hybrids_gbif_occ3 %>%
  filter(is.na(parent_scientificName)) %>%
  filter(!(scientificName %in% hybrids_w_parents$scientificName))
hybrids_gbif_occ3 <- bind_rows(hybrids_w_parents, hybrids_wo_parents)

hybrids_gbif_occ3 %<>%
  mutate(parent_scientificName = ifelse(
    is.na(parent_scientificName), parent_verbatimTaxonName, parent_scientificName
  )) %>%
  distinct()

all_hybrids8 <- all_hybrids7 %>%
  bind_rows(hybrids_gbif_occ3) %>%
  select(-verbatimTaxonName, -parent_verbatimTaxonName) %>%
  distinct() %>%
  filter(!grepl("× [a-zA-Z]$|× [a-aA-Z]{2}$", scientificName)) %>%
  filter(!grepl("^× [a-z]+$", scientificName))


all_hybrids8_clean <- all_hybrids8 %>%
  select(scientificName, parent_scientificName) %>%
  distinct()

all_hybrids8_clean_matched <- gbif.cmpfn(
  taxonName = unique(all_hybrids8_clean$scientificName),
  taxonID = seq_along(unique(all_hybrids8_clean$scientificName))
) %>%
  select(scientificName, usageKey, acceptedUsageKey, rank, status)

parent_hybrids8_clean_matched <- gbif.cmpfn(
  taxonName = unique(all_hybrids8_clean$parent_scientificName),
  taxonID = seq_along(unique(all_hybrids8_clean$parent_scientificName))
) %>%
  rename_with(~ paste0("parent_", .), .cols = everything()) %>%
  select(parent_scientificName, parent_usageKey, parent_acceptedUsageKey, parent_rank, parent_status)

all_hybrids8_clean %<>%
  mutate(parent_scientificName = ifelse(
    parent_scientificName == scientificName, NA, parent_scientificName
  )) %>%
  mutate(parent_scientificName = ifelse(
    gsub(" × ", " ", scientificName) == parent_scientificName, NA, parent_scientificName
  )) %>%
  left_join(all_hybrids8_clean_matched) %>%
  left_join(parent_hybrids8_clean_matched) %>%
  distinct()

hybrids_w_parents <- all_hybrids8_clean %>%
  filter(!is.na(parent_scientificName))
hybrids_wo_parents <- all_hybrids8_clean %>%
  filter(is.na(parent_scientificName)) %>%
  filter(!(scientificName %in% hybrids_w_parents$scientificName))
all_hybrids8_clean <- bind_rows(hybrids_w_parents, hybrids_wo_parents)

all_hybrids8_clean %<>% mutate(parent_scientificName = gsub("abelia abelia", NA, parent_scientificName))

all_hybrids8_clean %<>%
  mutate(scientificName = gsub("× [A-Z].* ", "× ", scientificName)) %>%
  distinct()
all_hybrids8_clean %<>% filter(scientificName != "Zygopetalum × f.")
all_hybrids8_clean %<>% filter(!grepl("^× [a-zA-Z.].*$", scientificName))
all_hybrids8_clean %<>% filter(!grepl("Miltoniopsis", scientificName))
all_hybrids8_clean %<>% filter(scientificName != "(Elymus trachycaulus × C)")
all_hybrids8_clean %<>% filter(scientificName != "(Agropyron × Alt")
all_hybrids8_clean %<>% filter(scientificName != "(Gentiana lawrencei × sp.")
all_hybrids8_clean %<>% mutate(scientificName = gsub("^\\(", "", scientificName))
all_hybrids8_clean %<>% mutate(scientificName = gsub("× ×", "×", scientificName))
all_hybrids8_clean %<>% distinct()
all_hybrids8_clean %>% write_csv2(here("hybrids_master_list.csv"))

all_hybrids8_clean %<>% mutate(acceptedUsageKey = ifelse(is.na(acceptedUsageKey), usageKey, acceptedUsageKey))
n_distinct(all_hybrids8_clean$acceptedUsageKey)

all_hybrids8_clean %>%
  filter(is.na(usageKey)) %>%
  pull(scientificName) %>%
  n_distinct()

all_hybrids8_clean %>%
  filter(is.na(parent_scientificName)) %>%
  pull(scientificName) %>%
  n_distinct()


# retrieve accepted names from GBIF -------------------------------------------------------------------

source(here("r_functions/fn_gbif_acceptedNames.R"))
gbif.cmpfn <- compiler::cmpfun(gbif.fn)

temp5 %<>%
  mutate(verbatimTaxonName = ifelse(is.na(verbatimTaxonName), scientificName, verbatimTaxonName)) %>%
  mutate(parent_verbatimTaxonName = ifelse(is.na(parent_verbatimTaxonName), parent_scientificName, parent_verbatimTaxonName)) %>%
  mutate(confidence = ifelse(matchType == "NONE", NA, confidence)) %>%
  mutate(matchType = ifelse(matchType == "NONE", NA, matchType)) %>%
  distinct()


gbif_matched <- gbif.cmpfn(
  taxonName = unique(temp5$scientificName),
  taxonID = seq_along(unique(temp5$scientificName))
)


# write_csv2(all_hybrids5, here("hybrids_ipni_and_floras.csv"))
# get chromosome numbers -----------------------------------------------------------------------
no_cores <- detectCores()
cl <- makeCluster(no_cores)
df <- temp3 %>%
  mutate(
    name_for_ccdb = str_squish(gsub("×", "", verbatimTaxonName)),
    parent_name_for_ccdb = str_squish(gsub("×", "", parent_verbatimTaxonName))
  )
ccdb_matches <- pblapply(unique(df$name_for_ccdb), chrom_counts, rank = "species", cl = cl)
ccdb_data <- ccdb_matches[sapply(ccdb_matches, nrow) > 0] %>%
  data.table::rbindlist()

parent_ccdb_matches <- pblapply(
  unique(df$parent_name_for_ccdb), chrom_counts,
  rank = "species", cl = cl
)
parents_ccdb_data <- parents_ccdb_matches[sapply(parents_ccdb_matches, nrow) > 0] %>%
  data.table::rbindlist()


# get POWO range info ------------------------------------------------------------------------
powo.fn <- function(id) {
  library(xml2)
  library(rvest)
  library(dplyr)
  myHTML <- try(read_html(paste0("http://powo.science.kew.org/taxon/", id)))
  suppressWarnings(if (class(myHTML) == "try-error") {
    range <- NA
  } else {
    range <- myHTML %>%
      html_elements(xpath = "//*[@id='distribution-listing']") %>%
      html_text() %>%
      stringr::str_squish()
    if (length(range) == 0) range <- NA else range <- range
  })
  df <- data.frame(id = id, range = range)
  return(df)
}

# ipni_id <- hybrids_ipni %>% pull(id) %>% unique()
# hybrid_range <- pblapply(ipni_id, powo.fn, cl = cl)
# hybrid_range_df <- hybrid_range %>%
#   bind_rows() %>%
#   filter(!is.na(range)) %>%
#   mutate(range = gsub("Doubtfully", "; Doubtfully", range)) %>%
#   mutate(range = gsub("Native", "; Native", range)) %>%
#   mutate(range = gsub("Introduced", "; Introduced", range)) %>%
#   mutate(range = gsub("Extinct", "; Extinct", range)) %>%
#   separate(range, into = letters[1:4], sep = ";") %>%
#   pivot_longer(cols = a:d, values_to = "range") %>%
#   filter(range != "" ) %>%
#   filter(!is.na(range)) %>%
#   mutate(range = str_trim(str_squish(range))) %>%
#   select(-name)

ipni_id <- hybrids_ipni %>%
  filter(parent_taxonID != "") %>%
  pull(parent_taxonID) %>%
  unique()
parent_range <- pblapply(ipni_id, powo.fn, cl = cl)
parent_range_df <- parent_range %>%
  bind_rows() %>%
  filter(!is.na(range)) %>%
  mutate(range = gsub("Doubtfully", "; Doubtfully", range)) %>%
  mutate(range = gsub("Native", "; Native", range)) %>%
  mutate(range = gsub("Introduced", "; Introduced", range)) %>%
  mutate(range = gsub("Extinct", "; Extinct", range)) %>%
  separate(range, into = letters[1:4], sep = ";") %>%
  pivot_longer(cols = a:d, values_to = "range") %>%
  filter(range != "") %>%
  filter(!is.na(range)) %>%
  mutate(range = str_trim(str_squish(range))) %>%
  select(-name) %>%
  rename(parent_taxonID = id) %>%
  filter(!grepl("Doubtfully|Extinct", range)) %>%
  separate(range, into = c("origin", "range"), sep = ":")

introduced_ranges <- parent_range_df %>%
  filter(grepl("Introduced", origin)) %>%
  rowwise() %>%
  mutate(region = str_split(range, ","), .keep = "unused") %>%
  unnest(region) %>%
  mutate(region = str_trim(region)) %>%
  mutate(origin = "Introduced")
# mutate(is_introduced = TRUE) %>%
# select(-origin) %>%
# left_join(hybrids_ipni) %>%
# select(id, is_introduced, region) %>%
# distinct()

native_ranges <- parent_range_df %>%
  filter(grepl("Native", origin)) %>%
  rowwise() %>%
  mutate(region = str_split(range, ","), .keep = "unused") %>%
  unnest(region) %>%
  mutate(region = str_trim(region)) %>%
  mutate(origin = "Native")
# mutate(is_native = TRUE) %>%
# select(-origin) %>%
# left_join(hybrids_ipni) %>%
# select(id, is_native, region) %>%
# distinct()

all_regions <- tibble(region = unique(c(native_ranges$region, introduced_ranges$region)))


ranges <- all_regions %>%
  left_join(native_ranges) %>%
  filter(!is.na(origin)) %>%
  bind_rows(
    all_regions %>%
      left_join(introduced_ranges) %>%
      filter(!is.na(origin))
  ) %>%
  inner_join(hybrids_ipni %>% select(id, parent_taxonID)) %>%
  distinct()

t <- ranges %>% count(origin, region, id)

ranges_list <- split(ranges, ranges$id)
ranges_list <- lapply(ranges_list, function(x) table(x$origin, x$region))
ranges_list2 <- sapply(ranges_list, function(x) length(table(colSums(x) > 1)) > 1)
ranges_list <- ranges_list[ranges_list2]
# df <- bind_rows(ranges_list)
df <- lapply(ranges_list, as.data.frame)
n <- names(ranges_list)

df2 <- bind_rows(df)
