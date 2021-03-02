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
  , nsim           = 2      ## Number of epidemic simulations for each filtering trajectory
  , loglik.thresh  = 2       ## Keep parameter sets with a likelihood within top X loglik units, to only fit with MLE, use 0
  , all.params     = FALSE
  , nparams        = 1
)

# model assessment ----
model_assess_params <- list( # general parameters for the model assessment simulations
     filtering.traj   = TRUE    
     , traj.trim_date = NA
     , sim_end_date   = as.Date("2020-07-01")
     , red_shelt.t    = 500 
     , ntraj          = 2
) 

model_asses_sims <- mdply(all_fits_trajs[5, ], 
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
source("intervention_sims.R")
# infected isolation sensitivity ---- 
source("inf_iso_sims.R")
# counterfactuals ----
source("counterfactual_sims.R")




