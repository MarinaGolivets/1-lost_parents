# Marina Golivets
# BiCIKL Hackathon (Sept 20-24, 2021, Meise)
# Extract hybrid and parent taxon names from floras


# load packages ------------------------------------------------------------------------------------

pkgs <- c("here", "tidyverse", "magrittr")
sapply(pkgs, require, character.only = TRUE)


# get hybrid names from FlorKart, German SL, Pladias, Swiss Flora ----------------------------------

hybrids_florkart <- read_csv2(here("data/florkart_taxa.csv")) %>%
  filter(isHybrid == TRUE) %>%
  filter(taxonRank == "species") %>%
  pull(taxonName, acceptedTaxonName)

hybrids_gsl <- readxl::read_xlsx(here("data/GermanSL1.5_replaceID.xlsx")) %>%
  janitor::clean_names() %>%
  filter(grepl("×|^x ", taxon_name, ignore.case = TRUE)) %>%
  filter(taxon_rank == "SPE") %>%
  mutate(name_author12 = gsub("-", "", name_author12)) %>%
  mutate(taxon_name = str_trim(paste(taxon_name, name_author12))) %>%
  pull(taxon_name, neuer_name)

hybrids_pladias <- read_csv2(here("data/pladias_taxa.csv")) %>%
  filter(isHybrid == TRUE) %>%
  filter(taxonRank == "species") %>%
  pull(taxonName, acceptedTaxonName)

hybrids_swissflora <- readxl::read_xlsx(
  here("data/Checklist_2017_simple_version_20181017.xlsx"),
  skip = 4
) %>%
  slice(2:n()) %>%
  filter(grepl("×|^x ", Taxonname, ignore.case = TRUE)) %>%
  pull(Taxonname)

hybrids_floras <- tibble(
  verbatim = c(hybrids_florkart, hybrids_gsl, hybrids_pladias, hybrids_swissflora) %>%
    unique()
) %>%
  mutate(edited = gsub("\\(purpurea × repens)", "purpurea × repens", verbatim)) %>%
  filter(!(edited %in% c("Salix 'Americana'", "Spiraea 'Arguta'"))) %>%
  mutate(edited = gsub("\\?", "", edited)) %>%
  mutate(edited = gsub(" x |×|^X ", " × ", edited)) %>%
  mutate(edited = str_trim(str_squish(edited))) %>%
  mutate(parent = str_split(edited, pattern = "×")) %>%
  unnest(cols = "parent") %>%
  mutate(parent = str_trim(str_squish(parent))) %>%
  rowwise() %>%
  mutate(genus = str_split(edited, pattern = " ", simplify = TRUE)[1]) %>%
  mutate(parent = case_when(
    grepl("^[a-z]", parent) ~ paste(genus, parent),
    grepl("^[A-Z]\\.", parent) ~
      paste(genus, gsub("^[A-Z]\\.", "", parent)),
    TRUE ~ parent
  )) %>%
  mutate(parent = ifelse( !grepl(" x |× ", edited, ignore.case = TRUE), NA, parent )) %>%
  mutate(parent = ifelse(grepl("× ", parent, ignore.case = TRUE), NA, parent)) %>%
  mutate(parent = ifelse(gsub(" × ", " ", edited) == parent, NA, parent)) %>%
  select(-genus) %>%
  filter(str_count(parent, pattern = " ") > 0 | is.na(parent)) %>%
  ungroup() %>%
  write_csv2(here("data/hybrids_floras.csv"))
