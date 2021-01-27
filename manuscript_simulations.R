all_fits <- c("output/Santa Clara_TRUE_FALSE_2020-04-22_2021-01-08_final.Rds")

# model assessment ----
model_assess_params <- 
  list(seed.val = 10001 
     , use.rds  = TRUE    
     , fit.E0   = TRUE    ## Was E0 also fit?
     # Whether to use filtering trajectories
     , filtering.traj = TRUE    ## simulate only from trajectories drawn from smoothing distribution
     , n.filter.traj  = 100      ## how many filtering trajectories to draw
     , nsim           = 100     ## Number of epidemic simulations for each filtering trajectory
     , n.particles    = 1000     ## how many particles to use in identifying the filtering trajectories
     # , sim_length     = 28      ## if using trajectories, how many additional days from the end of the data to run the simulation, if not using trajectories, how long total the simulation will be
     , sim_end_date   = as.Date("2020-07-01")
     , red_shelt.t    = 500     ## Time shelter in place changes to a reduced form (inf_iso or light) -- Now June 1 (day 76 of shelter in place orders)
     ## Intervention scenarios
     , inf_iso  = FALSE   ## Do we ever reduce from shelter in place to some form of strong/moderate social distancing?
     , light    = FALSE   ## Lightswitch method
     , loglik.thresh  = 2       ## Keep parameter sets with a likelihood within top X loglik units, to only fit with MLE, use 0
     , params.all = TRUE     ## Keep all fitted parameters above loglik thresh?...
) 

model_asses_sims <- alply(all_fits, 1, function(fit_name){
  source("./COVID_simulate_params.R", local = T)
  rds.name <- fit_name
  print(rds.name)
  with(model_assess_params, {
    # from rds name, identify the date of the last data used in it 
    last_data_date <- strsplit("output/Santa Clara_TRUE_FALSE_2020-04-22_2021-01-08_final.Rds", split = "_")[[1]][4] %>% as.Date()
    sim_length <- difftime(sim_end_date, last_data_date, units = "days") %>% as.numeric()
    source("./COVID_simulate.R", local = T)
    SEIR.sim.all %>% 
      mutate(last_data = last_data_date) %>% 
      return
  })
})

model_asses_sims[[1]] %>% 
  ggplot(aes(x = date, y = D_new, 
             group = interaction(paramset, mif2_iter, traj_set, .id),
             color = .id == "traj")) +
  facet_wrap(~paramset + mif2_iter) + 
  geom_line() + theme_bw()

# NPIs ----
NPI_general_params <- list(
  seed.val = 10001 
  , use.rds  = TRUE    
  , fit.E0   = TRUE    ## Was E0 also fit?
  # Whether to use filtering trajectories
  , filtering.traj = TRUE    ## simulate only from trajectories drawn from smoothing distribution
  , n.filter.traj  = 10      ## how many filtering trajectories to draw
  , nsim           = 10     ## Number of epidemic simulations for each filtering trajectory
  , n.particles    = 100     ## how many particles to use in identifying the filtering trajectories
  , sim_length     = 150      ## if using trajectories, how many additional days from the end of the data to run the simulation, if not using trajectories, how long total the simulation will be
  , loglik.thresh  = 2       ## Keep parameter sets with a likelihood within top X loglik units, to only fit with MLE, use 0
  , params.all = TRUE     ## Keep all fitted parameters above loglik thresh?...
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

# this is currently really inefficiency because we get the filtering trajectories separately for each set of intervention parameters even though they should be identical
intervention_sims <- ldply(NPI_intervention_params, function(int_params){
  source("./COVID_simulate_params.R", local = T)
  rds.name <- "output/Santa Clara_TRUE_FALSE_2020-04-22_2021-01-08_final.Rds"
  with(c(int_params, NPI_general_params), {
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
    with(c(int_params, NPI_general_params, 
           red_shelt.t          = 7
           , inf_iso            = TRUE
           , light              = FALSE), {
             test_and_isolate_s = inf_iso_level     
             test_and_isolate_m = inf_iso_level
             red_shelt.s        = background_soc_dist
             scenario_name = paste0("iso_", inf_iso_level, "soc_dist", background_soc_dist)
             print(scenario_name)
             source("./COVID_simulate.R", local = T)
             SEIR.sim.all %>% 
               mutate(scenario = scenario_name) %>% 
               return
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


