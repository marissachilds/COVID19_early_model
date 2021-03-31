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


NPI_general_params <- list(
  sim_length       = as.numeric(as.Date("2021-06-30") - as.Date("2020-04-22"))      
  , rds.name       = "output/Santa_Clara_TRUE_FALSE_2020-04-22_2021-02-05_final.Rds"
  , traj.file      = "output/Santa_Clara_TRUE_FALSE_2020-04-22_2021-02-05_final_filter_traj.Rds"
  , seed.val       = 10001 
  , nsim           = 25       ## Number of epidemic simulations for each filtering trajectory
  , ntraj          = 25
  , loglik.thresh  = 2        ## Keep parameter sets with a likelihood within top X loglik units, to only fit with MLE, use 0
)

NPI_intervention_params <- list(
  maintain       = list(red_shelt.t          = 1000     
                        , inf_iso            = FALSE   
                        , light              = FALSE
                        , scenario_name = "maintain")
  , test_isolate = list(red_shelt.t          = as.numeric(as.Date("2020-05-31") - as.Date("2020-04-22"))
                        , inf_iso            = TRUE
                        , light              = FALSE
                        , test_and_isolate_s = 0.2     
                        , test_and_isolate_m = 0.3    
                        , red_shelt.s        = 0.5
                        , scenario_name = "isolate")
  , lift         = list(red_shelt.t          = as.numeric(as.Date("2020-05-31") - as.Date("2020-04-22"))
                        , inf_iso            = FALSE   
                        , light              = FALSE
                        , red_shelt.s        = 0
                        , scenario_name = "lift")
  , light        = list(red_shelt.t          = as.numeric(as.Date("2020-05-31") - as.Date("2020-04-22"))
                        , inf_iso            = FALSE
                        , light              = TRUE
                        , red_shelt.crossed  = 0.20
                        , red_shelt.free     = 0.80  
                        , thresh_H.start     = 60      
                        , thresh_H.end       = 20
                        , scenario_name = "light")
)

intervention_traj <- ldply(NPI_intervention_params, function(int_params){
  print(int_params[["scenario_name"]])
  do.call(COVID_simulate, 
          args = c(int_params, NPI_general_params, 
                   nparams        = 1, 
                   last_date_only = FALSE)) %>%
    select(paramset, mif2_iter, traj_set, .id, D, H, H_new, I, I_new_sympt, day, date) %>% 
    mutate(scenario = int_params[["scenario_name"]]) %>%
    return
})

intervention_summary <- ldply(NPI_intervention_params[-4], function(int_params){
  print(int_params[["scenario_name"]])
  do.call(COVID_simulate, 
          args = c(int_params, NPI_general_params, 
                   last_date_only = FALSE, 
                   summary_fun  = function(sims_df){
                     sims_df %>%
                       filter(.id != "traj") %>%
                       rowwise() %>% 
                       mutate(inf_simul = I + H)  %>% 
                       group_by(paramset, mif2_iter, traj_set, .id) %>%
                       summarise(D_max = max(D), 
                                 D_new_max = max(D_new), 
                                 H_max = max(H),
                                 H_new_max = max(H_new),
                                 inf_max = max(inf_simul), 
                                 inf_max_date = date[which.max(inf_simul)],
				 I_new_sympt_max = max(I_new_sympt), 
				 I_new_sympt_max_date = date[which.max(I_new_sympt)],
                                 D_new_max_date = date[which.max(D_new)],
                                 .groups = "drop") %>% return}
                     )) %>%
    mutate(scenario = int_params[["scenario_name"]]) %>%
    return
})

saveRDS(list(summary = intervention_summary, 
             traj = intervention_traj),
        "output/intervention_sims.Rds")