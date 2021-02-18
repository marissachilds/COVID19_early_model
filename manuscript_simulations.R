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

# all fit names, trajectories names, and the associated fit dates
all_fits_trajs = data.frame(
  fit_name = c("output/Santa_Clara_TRUE_FALSE_2020-04-01_2021-02-08_final.Rds",
               "output/Santa_Clara_TRUE_FALSE_2020-04-08_2021-02-05_final.Rds",
               "output/Santa_Clara_TRUE_FALSE_2020-04-15_2021-02-05_final.Rds",
               "output/Santa_Clara_TRUE_FALSE_2020-04-22_2021-02-05_final.Rds",
               "output/Santa_Clara_TRUE_FALSE_2020-04-29_2021-02-08_final.Rds",
               "output/Santa_Clara_TRUE_FALSE_2020-05-06_2021-02-05_final.Rds",
               "output/Santa_Clara_TRUE_FALSE_2020-05-13_2021-02-05_final.Rds",
               "output/Santa_Clara_TRUE_FALSE_2020-05-20_2021-02-05_final.Rds",
               "output/Santa_Clara_TRUE_FALSE_2020-05-27_2021-02-05_final.Rds",
               "output/Santa_Clara_TRUE_FALSE_2020-06-03_2021-02-05_final.Rds",
               "output/Santa_Clara_TRUE_FALSE_2020-06-10_2021-02-05_final.Rds",
               "output/Santa_Clara_TRUE_FALSE_2020-06-17_2021-02-05_final.Rds",
               "output/Santa_Clara_TRUE_FALSE_2020-06-24_2021-02-05_final.Rds")) %>% 
  separate(fit_name, into = c(NA, NA, NA, NA, "fit_date", NA), sep = "_", remove = FALSE) %>% 
  mutate(fit_date = as.Date(fit_date),
         traj_name = gsub(".Rds", "_filter_traj.Rds", fit_name))

# some general parameters for simulating
sim_params <- list(
  seed.val = 10001 
  , nsim           = 10     ## Number of epidemic simulations for each filtering trajectory
  , loglik.thresh  = 2       ## Keep parameter sets with a likelihood within top X loglik units, to only fit with MLE, use 0
  , all.params = FALSE
  , nparams = 2
)

# model assessment ----
model_assess_params <- list( # general parameters for the model assessment simulations
     filtering.traj   = TRUE    
     , traj.trim_date = NA
     , sim_end_date   = as.Date("2020-07-01")
     , red_shelt.t    = 500 
     , ntraj = 2
) 

model_asses_sims <- mdply(all_fits_trajs[c(1,5,9), ], 
                          function(fit_name, fit_date, traj_name){
  print(fit_date)
  sim_length <- difftime(model_assess_params[["sim_end_date"]], fit_date, units = "days") %>% as.numeric()
  do.call(COVID_simulate, 
          args = c(model_assess_params, sim_params, 
                   sim_length = sim_length,
                   rds.name = fit_name, traj.file = traj_name)) %>% 
    mutate(fit_date = fit_date)  %>% return
})

model_asses_sims %>% 
  select(-ends_with("name")) %>% 
  ggplot(aes(x = date, y = deaths, # can't plot trajectories with deaths, only D_new 
             group = interaction(paramset, mif2_iter, traj_set, .id, fit_date))) +
  facet_wrap(~fit_date) + 
  geom_line(alpha = 0.2) + 
  geom_vline(aes(xintercept = fit_date)) +
  theme_bw()


# NPIs ----
NPI_general_params <- list(
  sim_length     = 300      ## if using trajectories, how many additional days from the end of the data to run the simulation, if not using trajectories, how long total the simulation will be
  , params.all     = FALSE     ## Keep all fitted parameters above loglik thresh?...
  , nparams        = 5  
)

NPI_intervention_params <- list(
  maintain       = list(red_shelt.t    = 500     
                        , inf_iso  = FALSE   
                        , light    = FALSE
                        , scenario_name = "maintain")
  , test_isolate = list(red_shelt.t          = 0
                        , inf_iso            = TRUE
                        , light              = FALSE
                        , test_and_isolate_s = 0.2     
                        , test_and_isolate_m = 0.2     
                        , red_shelt.s        = 0.5
                        , scenario_name = "isolate")
  , light        = list(red_shelt.t          = 0
                        , inf_iso            = FALSE
                        , light              = TRUE
                        , red_shelt.s        = 0.9
                        , thresh_H.start     = 100      
                        , thresh_H.end       = 50
                        , scenario_name = "light")
  , lift         = list(red_shelt.t = 0
                        , inf_iso  = FALSE   
                        , light    = FALSE
                        , red_shelt.s = 0
                        , scenario_name = "lift")
)

intervention_sims <- ldply(NPI_intervention_params, function(int_params){
  print(int_params[["scenario_name"]])
  do.call(COVID_simulate, 
          args = c(int_params, NPI_general_params, sim_params, 
                   rds.name  = "output/Santa_Clara_TRUE_FALSE_2020-04-29_2021-02-08_final.Rds",
                   traj.file = "output/Santa_Clara_TRUE_FALSE_2020-04-29_2021-02-08_final_filter_traj.Rds")) %>%
    mutate(scenario = int_params[["scenario_name"]]) %>%
    return
  })

intervention_sims %>% 
  filter(.id == 1) %>%
  ggplot(aes(x = date, y = H, 
             group = interaction(paramset, mif2_iter, 
                                 traj_set, .id, scenario),
             color = scenario)) +
  geom_line(alpha = 0.1) + 
  facet_wrap(~scenario, scales = "free") +
  scale_y_continuous(trans = "pseudo_log") + theme_bw()

intervention_sims %>% 
  filter(scenario == "light") %>% 
  filter(.id == 1) %>%
  select(paramset, mif2_iter, traj_set, .id, date, date, H, thresh_crossed) %>% 
  ggplot(aes(x = date, y = H, 
             color = thresh_crossed,
             group = interaction(paramset, mif2_iter, traj_set, .id))) + 
  geom_line(alpha = 0.1) 

# infected isolation sensitivity ---- 
inf_iso_sims <-
  expand.grid(inf_iso_level = c(0.1, 0.2, 0.3),
                            background_soc_dist = c(0.6, 0.5, 0.4)) %>% 
  mdply(function(inf_iso_level, background_soc_dist){

    scenario_name = paste0("iso_", inf_iso_level, "_soc_dist_", background_soc_dist)
    print(scenario_name)
    
    do.call(COVID_simulate, 
            args = c(NPI_general_params
                     , sim_params
                     , test_and_isolate_s = inf_iso_level     
                     , test_and_isolate_m = inf_iso_level
                     , red_shelt.t          = 7
                     , red_shelt.s        = background_soc_dist
                     , rds.name = "output/Santa_Clara_TRUE_FALSE_2020-04-29_2021-02-08_final.Rds"
                     , traj.file = "output/Santa_Clara_TRUE_FALSE_2020-04-29_2021-02-08_final_filter_traj.Rds"
                     , inf_iso            = TRUE
                     , light              = FALSE)) %>% 
               mutate(scenario = scenario_name) %>% 
               return
             # probably want to do summarizing in here to prevent this from getting too big
           })

inf_iso_sims %>% 
  filter(as.numeric(traj_set) <= 10) %>%
  ggplot(aes(x = date, y = I_new_sympt, 
             group = interaction(paramset, scenario, mif2_iter, traj_set, .id),
             color = scenario)) +
  facet_wrap(~scenario) + 
  geom_line() + 
  scale_y_continuous(trans = "pseudo_log") + theme_bw()

# counterfactuals ----
counterfactual_general_params <- list(
  params.all     = FALSE     ## Keep all fitted parameters above loglik thresh?...
  , nparams        = 5  
  , rds.name = "output/Santa_Clara_TRUE_FALSE_2020-04-29_2021-02-08_final.Rds"
  , nsim = 1 ## Number of epidemic simulations for each parameter set or filtering traj
  , traj.file = "output/Santa_Clara_TRUE_FALSE_2020-04-29_2021-02-08_final_filter_traj.Rds"
  , sim_length = 42      ## if using trajectories, how many additional days from the end of the data to run the simulation, if not using trajectories, how long total the simulation will be
  , traj.trim_date = as.Date("2020-03-16")
  , seed.val = 10001
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
                        , test_and_isolate_m = 0.2   
                        , test_and_isolate_s = 0.2   
                        , red_shelt.s        = 0.5
                        , scenario_name = "iso")
)


counterfactual_sims <- ldply(counterfactual_params, function(int_params){
  do.call(COVID_simulate, 
          args = c(counterfactual_general_params, int_params)) %>% 
    mutate(scenario = int_params[["scenario_name"]]) %>% return
})

counterfactual_sims %>% 
  ggplot(aes(x = date, y = D, 
             color = ifelse(.id == "traj", "traj", scenario),
             group = interaction(paramset, scenario, 
                                 mif2_iter, traj_set, .id))) + 
  geom_line(alpha = 0.5) + 
  geom_vline(xintercept = as.Date("2020-03-17")) +
  theme_bw()


counterfactual_sims %>%
  filter(date == max(date)) %>%
  select(traj_set, .id, date, paramset, mif2_iter, scenario, D) %>% 
  pivot_wider(names_from = scenario, values_from = D) %>% 
  mutate(delay_change = week_delay - observed, 
         iso_change = iso - observed) %>% 
  pivot_longer(ends_with("change")) %>% 
  ggplot(aes(x = value, group = name, 
             y = ..density..,
             color = name, 
             fill = name)) +
  geom_histogram(position = "identity", alpha = 0.5) + 
  geom_density(alpha = 0.5) + 
  geom_vline(xintercept = 0) +
  theme_bw()





