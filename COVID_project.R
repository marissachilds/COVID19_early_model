##################################################################
## Use pomp fits to simulate dynamics for various interventions ##
##################################################################

## Note: This script is set up to be able to be run from its own from a saved .Rds.
 ## Thus, packages are set to be loaded etc. (these lines do not need to be run if COVID_fit.R was just run)

####
## Parameters
####
set.seed(10001)
focal.county       <- "Santa Clara"
fitting            <- FALSE   ## Small change in pomp objects if fitting or simulating
## TRUE if COVID_fit previously run, FALSE if COVID_fit was just run and global environment is still full
use.rds            <- TRUE    
rds.name           <- "output/Santa Clara_TRUE_FALSE_2020-04-22_2021-01-08_final.Rds"
more.params.uncer  <- FALSE   ## Fit with more (FALSE) or fewer (TRUE) point estimates for a number of parameters
uncer.within.set   <- TRUE    ## Use parameter estimates from all mif2 runs for all parameter sets
nsim               <- 10     ## Number of epidemic simulations for each parameter set
fit.E0             <- TRUE    ## Was E0 also fit?
filtering.traj     <- TRUE    ## simulate only from trajectories drawn from smoothing distribution
n.filter.traj      <- 10     ## how many filtering trajectories to draw
n.particles        <- 500    ## how many particles to use in identifying the filtering trajectories
## Intervention scenarios. Only one at a time allowed right now! 
inf_iso            <- FALSE   ## Do we ever reduce from shelter in place to some form of strong/moderate social distancing?
test_and_isolate_s <- 0.2     ## Additional proportional reduction of severe cases under test and isolate
test_and_isolate_m <- 0.2     ## Additional proportional reduction of mild cases under test and isolate
light              <- FALSE   ## Lightswitch method
red_shelt.t        <- 76      ## Time shelter in place changes to a reduced form (inf_iso or light) -- Now June 1 (day 76 of shelter in place orders)
red_shelt.s        <- 0.5     ## New social dist strength after time red_shelt.t
thresh_H.start     <- 15      ## Threshold when lightswtich turns on (when we get higher than this)
thresh_H.end       <- 5       ## Threshold when lightswtich turns off (when we drop from above to this value)
# sim_length         <- 400     ## How many days to run the simulation
sim_length         <- 14      ## how many additional days from the end of the data to run the simulation
state.plot         <- "D"     ## State variable for plotting (Hospit [H], Death [D], or Cases [C])

loglik.thresh      <- 2       ## Keep parameter sets with a likelihood within top X loglik units, to only fit with MLE, use 0
params.all         <- TRUE   ## Keep all fitted parameters above loglik thresh?...
nparams            <- 200      ## ...if FALSE, pick the top X by loglik to use
plot.log10         <- FALSE   ## Plot on a log10 scale or not

## Required packages to run this code
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
lapply(needed_packages, require, character.only = TRUE)

## Theme for pretty plots
source("ggplot_theme.R")

## Load the previously saved fits
 ## If COVID_fit.R was just run, use parameers already stored in the global env 
if (use.rds) {
prev.fit           <- readRDS(rds.name)
variable_params    <- prev.fit[["variable_params"]] %>% 
  mutate(rowid = 1:n()) # add a row id to be able to match up with the data and covariate tables stored in data_covar
fixed_params       <- prev.fit[["fixed_params"]]
dynamics_summary   <- prev.fit[["dynamics_summary"]]
mif_traces         <- prev.fit[["mif_traces"]] 
data_covar         <- prev.fit[["data_covar"]]
pomp_components    <- prev.fit[["pomp_components"]]
}

## drop the rows that have 0s for likelihood (in case exited prematurely) 
## and keep only the best fits as defined by loglik
variable_params <- variable_params %>% 
  filter(log_lik != 0) %>% 
  filter(log_lik >= (max(log_lik) - loglik.thresh))

if (!params.all) {
  ## old option to pick a random X
  # variable_params <- variable_params[sample(1:nrow(variable_params), nparams), ]
  ## new option to subset to top X
  variable_params <- variable_params %>% arrange(desc(log_lik)) %>% slice(1:nparams)
}

####
## Simulate: loop over rows of variable_params
####

for (i in 1:nrow(variable_params)) {
print(i)
## define the parameters 
  full_params = c(fixed_params %>% 
                    filter(paramset == variable_params[[i, "paramset"]] & 
                             mif2_iter == variable_params[[i, "mif2_iter"]]) %>% 
                    select(param, value) %>% tibble::deframe(), 
                  variable_params[i,] %>% rename(beta0 = beta0est) %>% 
                    select(beta0, Ca, alpha, mu, delta, 
                           E_init, soc_dist_level_sip) %>% unlist)
  # need to set up new intervention table, for now, just using the old one
  # then set up pomp object for simulating, it will auto detect the starting states based on their names in the params vector (which will have _0) appended to the state name
  covid_simulate <- do.call(
    pomp, 
    # for the simulator, we don't want rinit
    args = c(pomp_components[names(pomp_components) != "rinit"],
             list(data = data_covar[[variable_params[i,]$rowid]]$county_data, 
                  covar = data_covar[[variable_params[i,]$rowid]]$covar_table)))
  
  # if using filtering trajectories, set up covid fitting object to pfilter
  covid_fitting <- do.call(pomp, 
                           args = c(pomp_components, 
                                    list(data = data_covar[[variable_params[i,]$rowid]]$county_data, 
                                         covar = data_covar[[variable_params[i,]$rowid]]$covar_table)))
  
  
  # pfilter to get filtering trajectories  
  trajs <- replicate(n.filter.traj, 
                     filter.traj(pfilter(covid_fitting, 
                                         params = full_params,
                                         filter.traj = TRUE, Np = n.particles))) 
  
  # for each filtering trajecory, simulate forward 
  SEIR.sim <- adply(.data = trajs[,1,last(covid_fitting@times),], 
                .margins = 2, 
                .fun = function(init_params){
                  names(init_params) = paste0(names(init_params), "_0")
                  simulate(covid_simulate,
                           nsim = nsim,
                           t0 = last(covid_fitting@times),
                           times = last(covid_fitting@times) + 1:sim_length,
                           format = "data.frame",
                           params = c(full_params, init_params))
                }) %>% 
    rename(traj_set = X1) %>% 
    mutate(paramset = variable_params[[i, "paramset"]],
           mif2_iter = variable_params[[i, "mif2_iter"]],
           date = round(as.Date(day, origin = variable_params[i, ]$sim_start)))
  
  # do we want to save the filtering trajectories too? 
     
## Stich together output
if (i == 1) {
 SEIR.sim.all <- SEIR.sim  
} else {
 SEIR.sim.all <- rbind(SEIR.sim.all, SEIR.sim) 
}

## Keep track of progress
if (((i / 20) %% 1) == 0) {
  print(paste(round(i / nrow(variable_params), 2)*100, "% Complete", sep = ""))
}

}

SEIR.sim.all %>%
  ggplot(aes(x = date, y = I, group = interaction(paramset, mif2_iter, traj_set, .id),
             color = interaction(paramset, mif2_iter))) +
  geom_line()

