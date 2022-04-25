############################################################
### Helper function to adjust exposure of the lgcp model ###
############################################################

adjust_exposure <- function(data, ips){
  
  # adjust exposure of the model dependent on tagging effort
  e_adj <- dat %>% as_tibble() %>%
    group_by(year_i, season) %>%
    summarise(n_ind = n_distinct(id),
              n_locs = n(),
              adj = n_ind/n_locs) %>%
    ungroup()
  
  # remember not all seasons and not all years....
  foo <- expand.grid(year_i = c(1, 2, 3, 5), season = c(1:4)) %>% tibble()
  foo <- foo %>% arrange(year_i, season) %>%
    left_join(e_adj) %>%
    arrange(year_i, season) %>%
    mutate(adj = replace_na(adj, 1)) # replace NA matches with 1 so weights aren't adjusted
  
  ips <- ips %>% 
    st_as_sf %>% 
    left_join(foo, by = c("year_i" = "year_i", "season" = "season")) %>%
    mutate(weight = weight * adj) %>%
    dplyr::select("season", "year_i", "weight") %>%
    as_Spatial()
  
  return(ips)
}
