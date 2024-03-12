wfo_h <- read_rds(here("hybrids-wfo-and-powo_0.rds"))



ipni_id <- wfo_h %>%
  select(wfo_id, contains("ipni_id_")) %>%
  unnest(c(ipni_id_ipni)) %>%
  unnest(c(ipni_id_ipni2)) %>%
  unnest(c(ipni_id_wfo)) %>%
  unnest(c(ipni_id_pow)) %>%
  pivot_longer(ipni_id_manual:ipni_id_pow, values_to = "ipni_id", 
               names_to = "source", values_drop_na = TRUE) %>%
  mutate(source = case_match(
    source,
    c("ipni_id_ipni", "ipni_id_ipni2") ~ "IPNI",
    "ipni_id_wfo" ~ "WFO",
    "ipni_id_pow" ~ "POWO",
    "ipni_id_manual" ~ "manual"
    )) %>%
  distinct() %>%
  arrange(source) %>%
  group_by(wfo_id, ipni_id) %>%
  summarise(source = str_c(source, collapse = "|"))

#---------------------------------------------------------------------------------------------------
cl <- makeCluster(detectCores() - 1)
doParallel::registerDoParallel(cl)

# query IPNI using IPNI ID as provided in WFO
ipni_ids <- unique(ipni_id$ipni_id) %>% .[!str_detect(., "-4$")]

q_id_ipni_l <- foreach(
  i = seq_along(ipni_ids),
  .packages = c("kewr"),
  .combine = "c",
  .errorhandling = "pass"
) %dopar% {
  out <- lookup_ipni(ipni_ids[i], type = "taxon")
  return(list(out))
}
write_rds(q_id_ipni_l, here("q2_id_ipni_l.rds"))

q_id_ipni_l <- read_rds(here("q2_id_ipni_l.rds"))

# check which IPNI IDs were skipped
q_id_ipni_l %>%
  .[names(.) != "call"] %>%
  .[sapply(., typeof) != "list"]

# unlist IPNI query results
parents1 <- q_id_ipni_l %>%
  .[sapply(., typeof) == "list"] %>%
  map(., ~ .[c("id", "hybridParents")] %>% compact()) %>%
  data.table::rbindlist(., fill = TRUE) %>%
  distinct() %>%
  rename(ipni_id = id) %>%
  mutate(parent12 = seq_along(ipni_id), .by = ipni_id) %>%
  unnest_longer(hybridParents) %>%
  filter(hybridParents_id %in% c("id", "wfoId", "name", "authors")) %>%
  pivot_wider(
    id_cols = c(ipni_id, parent12), 
    names_from = hybridParents_id, values_from = hybridParents)
  

parents2 <- q_id_ipni_l %>%
  .[sapply(., typeof) == "list"] %>%
  map(., ~ .[c("id", "remarks")]) %>%
  lapply(., function(x) x[sapply(x, length) != 0]) %>%
  data.table::rbindlist(., fill = TRUE) %>%
  filter(
    !is.na(remarks),
    str_detect(remarks, " X |hybr|cultiv")
    ) %>%
  distinct()

is_hybrid <- q_id_ipni_l %>%
  .[sapply(., typeof) == "list"] %>%
  map(., ~ .[c("id", "hybrid", "hybridGenus")]) %>%
  data.table::rbindlist(., fill = TRUE)
  
q_name_powo_l <- read_rds(here("q_name_powo_l.rds"))

# query POWO ---------------------------------------------------------------------------------------
powo_ids <- unique(ipni_id$ipni_id)

q_id_ipni_l <- foreach(
  i = seq_along(ipni_ids),
  .packages = c("kewr"),
  .combine = "c",
  .errorhandling = "pass"
) %dopar% {
  out <- lookup_powo(ipni_ids[i], distribution = TRUE)
  return(list(out))
}
