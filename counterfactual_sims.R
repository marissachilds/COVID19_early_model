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

counterfactual_general_params <- list(
  seed.val = 10001 
  , loglik.thresh  = 2       ## Keep parameter sets with a likelihood within top X loglik units, to only fit with MLE, use 0
  , rds.name = "output/fits/Santa_Clara_TRUE_FALSE_2020-04-22_2021-02-05_final.Rds"
  , nsim = 25 
  , ntraj = 25
  , traj.file = "output/filtering_trajectories/Santa_Clara_TRUE_FALSE_2020-04-22_2021-02-05_final_filter_traj.Rds"
  , sim_length = as.numeric(as.Date("2020-04-22") - as.Date("2020-03-16"))     
  , traj.trim_date = as.Date("2020-03-16")
  , counterfactual = TRUE
)

counterfactual_params <- list(
  obs       = list(red_shelt.t    = 500
                   , delay_days = 0
                   , scenario_name = "observed")
  , week_delay = list(red_shelt.t = 500     ## Time shelter in place changes to a reduced form (inf_iso or light), if using filtering trajectories, red_shelter.t is relative to the end of the data, otherwise its relative to when the shelter in place when into effect
                      , delay_days = 7
                      , scenario_name = "week_delay")
  , iso        = list(red_shelt.t = 0
                      , delay_days = 0
                      , inf_iso = TRUE
                      , test_and_isolate_m = 0.3   
                      , test_and_isolate_s = 0.2   
                      , red_shelt.s        = "sip"
                      , scenario_name = "iso")
)

counterfactual_sims_D <- ldply(counterfactual_params, function(int_params){
  do.call(COVID_simulate, 
          args = c(counterfactual_general_params, int_params,
                   last_date_only = TRUE)) %>% 
    select(paramset, mif2_iter, traj_set, .id, D, date) %>% 
    mutate(scenario = int_params[["scenario_name"]]) %>% return
})

counterfactual_sims_traj <- ldply(counterfactual_params, function(int_params){
  do.call(COVID_simulate, 
          args = c(counterfactual_general_params, int_params,
                   nparams        = 1,
                   last_date_only = FALSE)) %>% 
    select(paramset, mif2_iter, traj_set, .id, D,  I_new_sympt, day, date) %>% 
    mutate(scenario = int_params[["scenario_name"]]) %>% return
})


saveRDS(list(total_D = counterfactual_sims_D, 
             traj = counterfactual_sims_traj),
        "output/simulations/counterfactual_sims.Rds")