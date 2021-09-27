# Marina Golivets
# AlienScenarios WP5
# Function for standardizing taxon names against GBIF taxonomy


gbif.fn <- function(taxonName, taxonID) {
  library(dplyr)
  library(tidyr)

  no_cores <- parallel::detectCores()
  cl <- parallel::makeCluster(no_cores)
  all_matches <- pbapply::pblapply(
    taxonName, rgbif::name_backbone_verbose, kingdom = "plants", strict = TRUE, cl = cl
  )
  parallel::stopCluster(cl)

  altern <- lapply(
    all_matches,
    function(x) {
      y <- x$alternatives
      if (nrow(y) == 0) {
        y[1, 1] <- NA
        colnames(y) <- "usageKey"
      } else {
        y <- y
      }
      return(y)
    }
  ) %>%
    mapply(cbind, .,
      taxonName = taxonName, taxonID = taxonID,
      stringsAsFactors = FALSE, SIMPLIFY = FALSE
    ) %>%
    data.table::rbindlist(fill = TRUE) %>%
    filter(!is.na(usageKey)) %>%
    distinct()


  best_match <- lapply(all_matches, function(x) x$data) %>%
    mapply(cbind, .,
      taxonName = taxonName, taxonID = taxonID,
      stringsAsFactors = FALSE, SIMPLIFY = FALSE
    ) %>%
    data.table::rbindlist(fill = TRUE) %>%
    distinct()

  matched <- best_match %>% filter(!(matchType %in% c("NONE", "HIGHERRANK")))

  nonmatched <- best_match %>% filter(matchType %in% c("NONE", "HIGHERRANK"))

  matched_altern <- try(
    altern %>%
      filter(phylum == "Tracheophyta") %>%
      filter(confidence >= 0) %>%
      filter(taxonID %in% nonmatched$taxonID)
  )
  if (class(matched_altern)[1] == "try-error") {
    taxa_gbif <- matched %>%
      filter(rank != "GENUS")
  } else {
    taxa_gbif <- bind_rows(matched, matched_altern) %>%
      filter(rank != "GENUS")
  }

  accepted <- taxa_gbif %>%
    group_by(taxonID) %>%
    filter(status == "ACCEPTED") %>%
    filter(confidence == max(confidence)) %>%
    ungroup()

  synonyms <- taxa_gbif %>%
    group_by(taxonID) %>%
    summarise(hasAccept = length(unique(status == "ACCEPTED")) > 1) %>%
    full_join(taxa_gbif) %>%
    filter(hasAccept == FALSE) %>%
    filter(status == "SYNONYM") %>%
    group_by(taxonID) %>%
    filter(confidence == max(confidence)) %>%
    ungroup()

  doubtful <- taxa_gbif %>%
    group_by(taxonID) %>%
    summarise(hasAccept = length(unique(status == "ACCEPTED")) > 1) %>%
    full_join(taxa_gbif) %>%
    filter(hasAccept == FALSE) %>%
    group_by(taxonID) %>%
    filter(status == "DOUBTFUL") %>%
    filter(confidence == max(confidence)) %>%
    ungroup()

  taxa_gbif_edited <- bind_rows(accepted, synonyms, doubtful) %>%
    group_by(taxonID) %>%
    filter(confidence == max(confidence)) %>%
    filter(status != "NONE") %>%
    select(-hasAccept) %>%
    ungroup()

  return(taxa_gbif_edited)
}
