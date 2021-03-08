args<-commandArgs(TRUE)
#args<-"2020-04-01_2021-02-08"

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
  , "foreach"
  , "doParallel"
  , "stringr"
)


## load packages. Install all packages that return "FALSE"
lapply(needed_packages, library, character.only = TRUE)



# all fit names, trajectories names, and the associated fit dates
fit_name <- paste0("/home/users/marissac/COVID19_early_model/output/Santa_Clara_TRUE_FALSE_", args[1], "_final.Rds")
fit_date <- str_split(fit_name, "_")[[1]][7] %>%
	     as.Date()
traj_name <- gsub(".Rds", "_filter_traj.Rds", fit_name)         

source("COVID_simulate.R")         

# some general parameters for simulating
sim_params <- list(
  seed.val = 10001
  , nsim = 100     ## Number of epidemic simulations for each filtering$  
  , loglik.thresh  = 2       ## Keep parameter sets with a likelihood within top$  
  , all.params = TRUE
)

  
# model assessment ----
model_assess_params <- list( # general parameters for the model assessment simul$
  filtering.traj   = TRUE
  , traj.trim_date = NA
  , sim_end_date   = as.Date("2020-07-01")
  , red_shelt.t    = 500
  , ntraj = 100
)

sim_length <- difftime(model_assess_params[["sim_end_date"]], fit_date, units = "days") %>% as.numeric()

model_assess_sims <- COVID_simulate( filtering.traj   = TRUE
  , traj.trim_date = NA
  , red_shelt.t    = 500
  , nsim = 25     ## Number of epidemic simulations for each filtering$
  , loglik.thresh  = 2       ## Keep parameter sets with a likelihood within top$
  , all.params = TRUE
  , nparams=NA
  , sim_length = difftime(as.Date("2020-08-01"), fit_date, units = "days") %>% as.numeric()
  , ntraj=25,
               rds.name = fit_name, traj.file = traj_name) %>% 
               mutate(fit_date = fit_date) 
 

saveRDS(model_assess_sims, paste0("output/", fit_date, "_", Sys.Date(), "assess.Rds"))
