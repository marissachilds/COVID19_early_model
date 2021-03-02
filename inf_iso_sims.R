needed_packages <- c(
  "pomp"
  , "plyr"
  , "dplyr"
  , "ggplot2"
  , "magrittr"
  , "scales"
  , "lubridate"
  , "tidyr"
  , "data.table"
)

## load packages. Install all packages that return "FALSE"
lapply(needed_packages, library, character.only = TRUE)

source("COVID_simulate.R")

inf_iso_sims <-
  expand.grid(inf_iso_level = c(0.1, 0.2, 0.3),  # this is portion of remaining contacts for isolated infected 
              background_soc_dist = c(0.6, 0.5, 0.4, 0.3)) %>%  # this is scaling on beta, so 0.3 is most effective
  mdply(function(inf_iso_level, background_soc_dist){
    
    scenario_name = paste0("iso_", inf_iso_level, "_soc_dist_", background_soc_dist)
    print(scenario_name)
    
    do.call(COVID_simulate, 
            args = list(sim_length       = as.numeric(as.Date("2021-06-30") - as.Date("2020-04-22"))      
                        , seed.val       = 10001 
                        , nsim           = 25      
                        , ntraj          = 25
                        , loglik.thresh  = 2       
                        , test_and_isolate_s = inf_iso_level     
                        , test_and_isolate_m = inf_iso_level
                        , red_shelt.t        = as.numeric(as.Date("2020-05-31") - as.Date("2020-04-22")) # number of days that existing SIP should last
                        , red_shelt.s        = background_soc_dist
                        , rds.name  = "output/Santa_Clara_TRUE_FALSE_2020-04-22_2021-02-05_final.Rds"
                        , traj.file = "output/Santa_Clara_TRUE_FALSE_2020-04-22_2021-02-05_final_filter_traj.Rds"
                        , inf_iso            = TRUE
                        , light              = FALSE
                        , last_date_only     = TRUE)) %>% 
      select(paramset, mif2_iter, traj_set, .id, D, date) %>% 
      mutate(scenario = scenario_name) %>% 
      return
  })

inf_iso_sims_comp <-
  do.call(COVID_simulate, 
          args = list(sim_length       = as.numeric(as.Date("2021-06-30") - as.Date("2020-04-22"))      
                      , seed.val       = 10001 
                      , nsim           = 25      
                      , ntraj          = 25
                      , loglik.thresh  = 2       
                      , red_shelt.t        = 1000
                      , red_shelt.s        = "sip"
                      , rds.name  = "output/Santa_Clara_TRUE_FALSE_2020-04-22_2021-02-05_final.Rds"
                      , traj.file = "output/Santa_Clara_TRUE_FALSE_2020-04-22_2021-02-05_final_filter_traj.Rds"
                      , inf_iso            = FALSE
                      , light              = FALSE
                      , last_date_only     = TRUE)) %>% 
  select(paramset, mif2_iter, traj_set, .id, D, date) %>% 
  mutate(scenario = "maintain") 

inf_iso_sims %<>% rbind(inf_iso_sims_comp %>% 
                          mutate(inf_iso_level = NA, 
                                 background_soc_dist = NA))

saveRDS(inf_iso_sims, "output/inf_iso_sims_D.Rds")