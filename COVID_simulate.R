##################################################################
## Use pomp fits to simulate dynamics for various interventions ##
##################################################################

## Note: This script is set up to be able to be run from its own from a saved .Rds.
 ## Thus, packages are set to be loaded etc. (these lines do not need to be run if COVID_fit.R was just run)

# the script can either simulate forward from the end of the data using the 
# filtering trajectories or simulate forward from t0 

# the parameters that must be set for simulation are in the COVID_simulate_params.R file

library(doParallel)

COVID_simulate <- function(rds.name, 
                           nsim = 1, ## Number of epidemic simulations for each parameter set or filtering traj
                           traj.file = NULL, # NULL results in not using filtering trajectories
                           ntraj = NA, # if not NA, number of filtering trajectories to use
                           sim_length = 28,      ## if using trajectories, how many additional days from the end of the data to run the simulation, if not using trajectories, how long total the simulation will be
                           traj.trim_date = NA,
                           red_shelt.t,     ## Time shelter in place changes to a reduced form (inf_iso or light), if using filtering trajectories, red_shelter.t is relative to the end of the data, otherwise its relative to when the shelter in place when into effect
                           inf_iso = FALSE,   ## Do we ever go to infected isolation?
                           test_and_isolate_s = NA,     ## Additional proportional reduction of severe cases under test and isolate
                           test_and_isolate_m = NA,     ## Additional proportional reduction of mild cases under test and isolate
                           light = FALSE, # Do we ever do to lightswitch? 
                           red_shelt.s = NA,     ## New social dist strength after time red_shelt.t
                           red_shelt.crossed = NA, ## only used for lightswitch, distancing when hospitals are above a threshold
                           red_shelt.free    = NA, ## only used for lightswitch, distancing when hospitals are below a threshold
                           thresh_H.start = NA, ## Threshold when lightswtich turns on (when we get higher than this)
                           thresh_H.end = NA, ## Threshold when lightswtich turns on (when we get higher than this)
                           loglik.thresh = 2,       ## Keep parameter sets with a likelihood within top X loglik units, to only fit with MLE, use 0
                           nparams = NA,  # if no NA, the number of parameter sets to keep
                           seed.val = NULL, 
                           counterfactual = FALSE,
                           delay_days = NA,
                           last_date_only = FALSE,
                           ...){     ## ...if FALSE, pick the top X by loglik to use
                           
if(!is.null(seed.val)) {set.seed(seed.val)}

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

## Load the previously saved fits
 ## If COVID_fit.R was just run, use parameers already stored in the global env 
prev.fit           <- readRDS(rds.name)
variable_params    <- prev.fit[["variable_params"]] %>% 
  mutate(rowid = 1:n()) # add a row id to be able to match up with the data and covariate tables stored in data_covar and the trajectories if they get read in
fixed_params       <- prev.fit[["fixed_params"]]
dynamics_summary   <- prev.fit[["dynamics_summary"]]
mif_traces         <- prev.fit[["mif_traces"]] 
data_covar         <- prev.fit[["data_covar"]]
pomp_components    <- prev.fit[["pomp_components"]]

if(!is.null(traj.file)){
  trajs <- readRDS(traj.file)
}

## drop the rows that have 0s for likelihood (in case exited prematurely) 
## and keep only the best fits as defined by loglik
variable_params <- variable_params %>% 
  filter(log_lik != 0) %>% 
  filter(log_lik >= (max(log_lik) - loglik.thresh))

if (!is.na(nparams)) {
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
  
  if(!is.null(traj.file)){
    
    current_traj <- trajs[[variable_params[i,]$rowid]]
    traj_times <- as.numeric(dimnames(current_traj)[[2]])
    
    if(!is.na(ntraj)){ current_traj <- current_traj[,,1:min(dim(current_traj)[3], ntraj)]}
    
    # if doing a counterfactual, trim the trajectories to a certain date (if that date falls before the end of the trajectories)
    if(!is.na(traj.trim_date)){
      traj.trim_index <- which(traj_times == (data_covar[[variable_params[i,]$rowid]]$county_data %>% filter(date == traj.trim_date) %>% pull(day) %>% first))
      if(length(traj.trim_index) > 0){
        traj_times <- traj_times[1:traj.trim_index]
        current_traj <- current_traj[,1:traj.trim_index,]
      }}
    
    # set up an intervention table for simulating after the data end, using the previous one up until the data ends
    intervention.forecast = data_covar[[variable_params[i,]$rowid]]$covar_table
    # intervention table that will get added onto the previous one at the point where the data end
    intervention.forecast_cont <- with(variable_params[i,], {
      red_shelt.t = min(sim_length, red_shelt.t) # limit the SIP intervention to the length of the simulation
      covariate_table(day = last(traj_times) + 1:sim_length
                      ,intervention = c(
                        # continut SIP for red_shelt.t days
                        {if(counterfactual){c(rep(1, delay_days), rep(2, red_shelt.t - delay_days))} else {rep(2, red_shelt.t)}} 
                        # Post SIP close
                        , {if (!inf_iso & !light) { # no intervention
                          rep(0, max(sim_length - red_shelt.t, 0))
                        } else if (inf_iso & !light) { # infected isolation with background social
                          rep(1, max(sim_length - red_shelt.t, 0))
                        } else if (!inf_iso & light) { # lightswitch
                          rep(3, max(sim_length - red_shelt.t, 0))
                        }}
                      )
                      , thresh_int_level = if (!light) {
                        rep(soc_dist_level_sip, sim_length)
                      } else {
                        rep(red_shelt.crossed, sim_length)
                      }
                      , back_int_level   = if (!light) {
                        rep(red_shelt.s, sim_length)
                      } else {
                        rep(red_shelt.free, sim_length) 
                      }
                      , isolation = {if (!inf_iso) { rep(0, sim_length)
                                    } else { c(rep(0, red_shelt.t), rep(1,  max(sim_length - red_shelt.t, 0)))}}
                      , iso_severe_level = rep(test_and_isolate_s, sim_length)      # % of contats that severe cases maintain
                      , iso_mild_level   = rep(test_and_isolate_m, sim_length)      # % of contats that mild cases maintain
                      , soc_dist_level_wfh = {
                        if (!inf_iso & !light) {
                          rep(soc_dist_level_wfh, sim_length)
                        } else { # if there is infected isolation, go to the reduced social distancing level when it goes into effect
                          c(
                            rep(soc_dist_level_wfh, red_shelt.t)
                            , rep(red_shelt.s, max(sim_length - red_shelt.t, 0))
                          )
                        }
                      }
                      , thresh_H_start   = rep(thresh_H.start, sim_length) 
                      , thresh_H_end     = rep(thresh_H.end, sim_length)
                      , order            = "constant"
                      , times            = "day"  
                      )
      })
    # force new intervention table rowname order to match those of the existing intervention table
    intervention.forecast_cont@table <- intervention.forecast_cont@table[match(rownames(intervention.forecast@table), rownames(intervention.forecast_cont@table)), ]
    # replace the end of the old intervention table with this new one 
    cut_point = which(intervention.forecast@times == last(traj_times))
    intervention.forecast@times <- c(intervention.forecast@times[1:cut_point], intervention.forecast_cont@times)
    intervention.forecast@table <- cbind(intervention.forecast@table[,1:cut_point], intervention.forecast_cont@table)
    
  } else{
    
    # if just simulating, set up the intervention table without concern for the existing one
    # set up a new intervention table
    intervention.forecast <- with(variable_params[i, ], {
      covariate_table(
        day              = 1:sim_length
        , intervention     = c(
          # No intervention until intervention start time
          rep(0, int_start1 - sim_start)                   
          # Intervention style 1
          , rep(1, int_length1)
          # Intervention style 2
          , rep(2, red_shelt.t)
          # Post intervention close
          , {
            if (!inf_iso & !light) {
              rep(0, sim_length - (int_start2 - sim_start) - red_shelt.t)
            } else if (inf_iso & !light) {
              rep(1, sim_length - (int_start2 - sim_start) - red_shelt.t)  
            } else if (!inf_iso & light) {
              rep(3, sim_length - (int_start2 - sim_start) - red_shelt.t)        
            }
          }
        )
        , thresh_int_level = rep(soc_dist_level_sip, sim_length)
        , back_int_level   = rep(red_shelt.s, sim_length)
        , isolation = { 
          if (!inf_iso) { rep(0, sim_length)
          } else { c(rep(0, int_start2 + red_shelt.t - sim_start)  
                     , rep(1, sim_length - (int_start2 + red_shelt.t - sim_start)))
          }}
        , iso_severe_level = rep(test_and_isolate_s, sim_length)      # % of contats that severe cases maintain
        , iso_mild_level   = rep(test_and_isolate_m, sim_length)   # % of contats that mild cases maintain
        , soc_dist_level_wfh = {
          if (!inf_iso) {
            rep(soc_dist_level_sip, sim_length)
          } else {
            c(
              rep(soc_dist_level_sip, (int_start1 - sim_start) + int_length1 + red_shelt.t)
              , rep(red_shelt.s, sim_length - ((int_start1 - sim_start) + int_length1 + red_shelt.t))
            )
          }
        }
        , thresh_H_start   = rep(thresh_H.start, sim_length) 
        , thresh_H_end     = rep(thresh_H.end, sim_length)
        , order            = "constant"
        , times            = "day"
      )
      
    })
    
  }
  
  # then set up pomp object for simulating, 
  # for filtering trajectories, it will auto detect the starting states based on their names in the params vector (which will have _0) appended to the state name
  covid_simulate <- do.call(
    pomp, 
    # for the simulator, we don't want rinit
    args = c({if(!is.null(traj.file)) pomp_components[names(pomp_components) != "rinit"] else pomp_components},
             list(data = data_covar[[variable_params[i,]$rowid]]$county_data, 
                  covar = intervention.forecast)
             ))
  
  # simulate, either forward from the filtering trajectory, or from t0
  if(!is.null(traj.file)){
    SEIR.sim <- adply(.data = 1:dim(current_traj)[3], 
                      .margins = 1, 
                      .fun = function(j){
                        init_params = current_traj[,length(traj_times), j]
                        names(init_params) = paste0(names(init_params), "_0")
                        simulate(covid_simulate,
                                 nsim = nsim,
                                 seed = seed.val, 
                                 t0 = last(traj_times),
                                 times = last(traj_times) + 1:sim_length,
                                 format = "data.frame",
                                 params = c(full_params, init_params)) %>% 
                          # add on the filtering trajectories too
                          rbind(current_traj[,,j] %>% t %>% 
                                  as.data.frame() %>% 
                                  cbind(day = traj_times,
                                        .id = "traj", 
                                        deaths = NA, date = NA)) %>% 
                          return
                      }) %>% 
      rename(traj_set = X1) 
    
  } else{
    SEIR.sim <- simulate(covid_simulate,
                         nsim = nsim,
                         seed = seed.val, 
                         times = intervention.forecast@times,
                         format = "data.frame",
                         params = full_params)
  }
  
  
 SEIR.sim %<>% 
   mutate(paramset = variable_params[[i, "paramset"]],
          mif2_iter = variable_params[[i, "mif2_iter"]],
          date = round(as.Date(day, origin = variable_params[i, ]$sim_start)))
 
if (last_date_only) {
 SEIR.sim %<>% filter(date == max(date))
}
     
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

return(SEIR.sim.all)

}
