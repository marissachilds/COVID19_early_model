fits_trajs <- list(list(rds.name = "output/Santa_Clara_TRUE_FALSE_2020-04-01_2021-02-08_final.Rds",
                        traj.file = "output/Santa_Clara_TRUE_FALSE_2020-04-01_2021-02-08_final_filter_traj.Rds"),
                   list(rds.name = "output/Santa_Clara_TRUE_FALSE_2020-04-08_2021-02-05_final.Rds",
                        traj.file = "output/Santa_Clara_TRUE_FALSE_2020-04-08_2021-02-05_final_filter_traj.Rds"))

sim_params <- list(
  seed.val = 10001 
  , use.rds  = TRUE    
  , fit.E0   = TRUE    ## Was E0 also fit?
  , nsim           = 1     ## Number of epidemic simulations for each filtering trajectory
  , loglik.thresh  = 2       ## Keep parameter sets with a likelihood within top X loglik units, to only fit with MLE, use 0
  , params.all = TRUE     ## Keep all fitted parameters above loglik thresh?...
)

# model assessment ----
model_assess_params <- 
  list(
     filtering.traj   = TRUE    ## simulate only from trajectories drawn from smoothing distribution
     , traj.trim_date = NA
     , sim_end_date   = as.Date("2020-05-01")
     , red_shelt.t    = 500   
     , inf_iso        = FALSE   ## Do we ever reduce from shelter in place to some form of strong/moderate social distancing?
     , light          = FALSE   ## Lightswitch method
) 

model_asses_sims <- adply(fits_trajs, 1, function(file_names){
  source("./COVID_simulate_params.R", local = T)
  with(c(model_assess_params, sim_params, file_names), {
    print(rds.name)
    # from rds name, identify the date of the last data used in it
    last_data_date <- strsplit(rds.name, split = "_")[[1]][5] %>% as.Date()
    sim_length <- difftime(sim_end_date, last_data_date, units = "days") %>% as.numeric()
    source("./COVID_simulate.R", local = T)
    SEIR.sim.all %>%
      mutate(last_data = last_data_date) %>%
      return
  })
})

model_asses_sims %>% 
  ggplot(aes(x = date, y = D_new, 
             group = interaction(paramset, mif2_iter, traj_set, .id),
             color = interaction(paramset, mif2_iter))) +
  facet_wrap(~last_data) + 
  geom_line() + 
  geom_vline(aes(xintercept = last_data)) +
  theme_bw()

# NPIs ----
NPI_general_params <- list(
  filtering.traj = TRUE    ## simulate only from trajectories drawn from smoothing distribution
  , sim_length     = 150      ## if using trajectories, how many additional days from the end of the data to run the simulation, if not using trajectories, how long total the simulation will be
)

NPI_intervention_params <- list(
  maintain = list(red_shelt.t    = 300     
                  , inf_iso  = FALSE   
                  , light    = FALSE
                  , scenario_name = "maintain")
  , test_isolate = list(red_shelt.t          = 7
                        , inf_iso            = TRUE
                        , light              = FALSE
                        , test_and_isolate_s = 0.2     
                        , test_and_isolate_m = 0.2     
                        , red_shelt.s        = 0.5
                        , scenario_name = "isolate")
  , light = list(red_shelt.t          = 7
                , inf_iso            = FALSE
                , light              = TRUE
                , red_shelt.s        = 0.5 # social distancing level when light switch is off
                , thresh_H.start     = 600      
                , thresh_H.end       = 300
                , scenario_name = "light")
  , lift = list(red_shelt.t = 7
                , inf_iso  = FALSE   
                , light    = FALSE
                , red_shelt.s = 0
                , scenario_name = "lift")
)

# this is currently really inefficient because we get the filtering trajectories separately for each set of intervention parameters even though they should be identical
intervention_sims <- ldply(NPI_intervention_params, function(int_params){
  source("./COVID_simulate_params.R", local = T)
  rds.name  <- "output/Santa Clara_TRUE_FALSE_2020-06-24_2021-01-26_final.Rds"
  traj.file <- "output/Santa Clara_TRUE_FALSE_2020-06-24_2021-01-26_final_filter_traj.Rds"
  with(c(int_params, NPI_general_params, sim_params), {
    print(scenario_name)
    source("./COVID_simulate.R", local = T)
    SEIR.sim.all %>% 
      mutate(scenario = scenario_name) %>% 
      return
  })
})

intervention_sims %>% 
  # filter(scenario != "lift") %>% 
  ggplot(aes(x = date, y = H, 
             group = interaction(paramset, mif2_iter, traj_set, .id),
             color = scenario)) +
  facet_wrap(~scenario) + 
  geom_line() + 
  scale_y_continuous(trans = "pseudo_log") + theme_bw()

# infected isolation sensetivity --- 
inf_iso_sims <-
  expand.grid(inf_iso_level = c(0.1, 0.2, 0.3),
                            background_soc_dist = c(0.6, 0.5, 0.4)) %>% 
    {split(., f = 1:nrow(.))} %>% 
    lapply(as.list) %>% 
  ldply(function(int_params){
    source("./COVID_simulate_params.R", local = T)
    rds.name <- "output/Santa Clara_TRUE_FALSE_2020-04-22_2021-01-08_final.Rds"
    traj.file <- "output/Santa Clara_TRUE_FALSE_2020-06-24_2021-01-26_final_filter_traj.Rds"
    with(c(int_params, NPI_general_params, 
           sim_params,
           red_shelt.t          = 7
           , inf_iso            = TRUE
           , light              = FALSE), {
             test_and_isolate_s = inf_iso_level     
             test_and_isolate_m = inf_iso_level
             red_shelt.s        = background_soc_dist
             scenario_name = paste0("iso_", inf_iso_level, "_soc_dist_", background_soc_dist)
             print(scenario_name)
             source("./COVID_simulate.R", local = T)
             SEIR.sim.all %>% 
               mutate(scenario = scenario_name) %>% 
               return
             # probably want to do summarizing in here to prevent this from getting too big
           })
  })

inf_iso_sims %>% 
  ggplot(aes(x = date, y = H, 
             group = interaction(paramset, mif2_iter, traj_set, .id),
             color = scenario)) +
  facet_wrap(~scenario) + 
  geom_line() + 
  scale_y_continuous(trans = "pseudo_log") + theme_bw()

# counterfactuals ----


