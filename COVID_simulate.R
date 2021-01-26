##################################################################
## Use pomp fits to simulate dynamics for various interventions ##
##################################################################

## Note: This script is set up to be able to be run from its own from a saved .Rds.
 ## Thus, packages are set to be loaded etc. (these lines do not need to be run if COVID_fit.R was just run)

# the script can either simulate forward from the end of the data using the 
# filtering trajectories, in which case (currnetly) no interventions can be 
# applied, or simulate forward from t0 in which case multiple intervention scenarios are possible 

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
nsim               <- 10      ## Number of epidemic simulations for each parameter set
fit.E0             <- TRUE    ## Was E0 also fit?

# Whether to use filtering trajectories
filtering.traj     <- TRUE    ## simulate only from trajectories drawn from smoothing distribution
n.filter.traj      <- 10      ## how many filtering trajectories to draw
n.particles        <- 500     ## how many particles to use in identifying the filtering trajectories
sim_length         <- 28      ## if using trajectories, how many additional days from the end of the data to run the simulation, if not using trajectories, how long total the simulation will be

## Intervention scenarios. Only one at a time allowed right now! 
## if using filtering trajectories, red_shelter.t is relative to the end of the data, otherwise its relative to when the shelter in place when into effect
inf_iso            <- FALSE   ## Do we ever reduce from shelter in place to some form of strong/moderate social distancing?
test_and_isolate_s <- 0.2     ## Additional proportional reduction of severe cases under test and isolate
test_and_isolate_m <- 0.2     ## Additional proportional reduction of mild cases under test and isolate
red_shelt.t        <- 500     ## Time shelter in place changes to a reduced form (inf_iso or light) -- Now June 1 (day 76 of shelter in place orders)
red_shelt.s        <- 0.5     ## New social dist strength after time red_shelt.t
light              <- FALSE   ## Lightswitch method
thresh_H.start     <- 15      ## Threshold when lightswtich turns on (when we get higher than this)
thresh_H.end       <- 5       ## Threshold when lightswtich turns off (when we drop from above to this value)
# sim_length         <- 400     ## How many days to run the simulation
state.plot         <- "D"     ## State variable for plotting (Hospit [H], Death [D], or Cases [C])

loglik.thresh      <- 10       ## Keep parameter sets with a likelihood within top X loglik units, to only fit with MLE, use 0
params.all         <- TRUE     ## Keep all fitted parameters above loglik thresh?...
nparams            <- 200      ## ...if FALSE, pick the top X by loglik to use
plot.log10         <- FALSE    ## Plot on a log10 scale or not

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
# source("ggplot_theme.R")

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

## Adjust variable params for the simulation scenario
variable_params <- variable_params %>%
  mutate(   
    int_length2        = red_shelt.t
    , iso_start          = int_start2 + red_shelt.t
    , soc_dist_level_red = red_shelt.s
  )


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
  
  if(filtering.traj){
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
    
    # set up an intervention table for simulating after the data end, using the previous one up until the data ends
    intervention.forecast = data_covar[[variable_params[i,]$rowid]]$covar_table
    # intervention table that will get added onto the previous one at the point where the data end
    intervention.forecast_cont <- with(variable_params[i,], {
      int_length2 = min(sim_length, int_length2) # limit the SIP intervention to the length of the simulation
      covariate_table(day = last(covid_fitting@times) + 1:sim_length
                      ,intervention = c(
                        # continut SIP for red_shelt.t days
                        rep(2, int_length2)
                        # Post SIP close
                        , {if (!inf_iso & !light) { # no intervention
                          rep(0, max(sim_length - int_length2, 0))
                        } else if (inf_iso & !light) { # infected isolation
                          rep(1, max(sim_length - int_length2, 0))
                        } else if (!inf_iso & light) { # lightswitch
                          rep(3, max(sim_length - int_length2, 0))
                        }}
                      )
                      , thresh_int_level = rep(soc_dist_level_sip, sim_length)
                      , back_int_level   = rep(red_shelt.s, sim_length)
                      , isolation = {if (!inf_iso) { rep(0, sim_length)
                                    } else { c(rep(0, int_length2), rep(1,  max(sim_length - int_length2, 0)))
                        }}
                      , iso_severe_level = rep(test_and_isolate_s, sim_length)      # % of contats that severe cases maintain
                      , iso_mild_level   = rep(test_and_isolate_m, sim_length)   # % of contats that mild cases maintain
                      , soc_dist_level_wfh = {
                        if (!inf_iso) {
                          rep(soc_dist_level_wfh, sim_length)
                        } else { # if there is infected isolation, go to the reduced social distancing level when it goes into effect
                          c(
                            rep(soc_dist_level_wfh, int_length2)
                            , rep(soc_dist_level_red, max(sim_length - int_length2, 0))
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
    cut_point = which(intervention.forecast@times == last(covid_fitting@times))
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
          , rep(2, int_length2)
          # Post intervention close
          , {
            if (!inf_iso & !light) {
              rep(0, sim_length - (int_start2 - sim_start) - int_length2)
            } else if (inf_iso & !light) {
              rep(2, sim_length - (int_start2 - sim_start) - int_length2)  
            } else if (!inf_iso & light) {
              rep(3, sim_length - (int_start2 - sim_start) - int_length2)        
            }
          }
        )
        , thresh_int_level = rep(soc_dist_level_sip, sim_length)
        , back_int_level   = rep(red_shelt.s, sim_length)
        , isolation = { 
          if (!inf_iso) { rep(0, sim_length)
          } else { c(rep(0, iso_start - sim_start)  
                     , rep(1, sim_length - (iso_start - sim_start)))
          }}
        , iso_severe_level = rep(test_and_isolate_s, sim_length)      # % of contats that severe cases maintain
        , iso_mild_level   = rep(test_and_isolate_m, sim_length)   # % of contats that mild cases maintain
        , soc_dist_level_wfh = {
          if (!inf_iso) {
            rep(soc_dist_level_sip, sim_length)
          } else {
            c(
              rep(soc_dist_level_sip, (int_start1 - sim_start) + int_length1 + int_length2)
              , rep(soc_dist_level_red, sim_length - ((int_start1 - sim_start) + int_length1 + int_length2))
            )
          }
        }
        # , soc_dist_level_wfh = rep(soc_dist_level_wfh, sim_length) 
        # , soc_dist_level_sip = {
        #   if (!inf_iso) {
        #     rep(soc_dist_level_sip, sim_length) 
        #   } else {
        #     c(
        #       rep(soc_dist_level_sip, (int_start1 - sim_start) + int_length1 + int_length2)
        #       , rep(soc_dist_level_red, sim_length - ((int_start1 - sim_start) + int_length1 + int_length2))
        #     )
        #   }
        # }
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
    args = c({if(filtering.traj) pomp_components[names(pomp_components) != "rinit"] else pomp_components},
             list(data = data_covar[[variable_params[i,]$rowid]]$county_data, 
                  covar = intervention.forecast)
             ))
  
  # for each filtering trajecory, simulate forward 
  if(filtering.traj){
    SEIR.sim <- adply(.data = 1:dim(trajs)[4], 
                      .margins = 1, 
                      .fun = function(j){
                        init_params = trajs[,1,last(covid_fitting@times), j]
                        names(init_params) = paste0(names(init_params), "_0")
                        simulate(covid_simulate,
                                 nsim = nsim,
                                 t0 = last(covid_fitting@times),
                                 times = last(covid_fitting@times) + 1:sim_length,
                                 format = "data.frame",
                                 params = c(full_params, init_params)) %>% 
                          # add on the filtering trajectories too
                          rbind(trajs[,1,, j] %>% t %>% 
                                  as.data.frame() %>% 
                                  cbind(day = covid_fitting@t0:last(covid_fitting@times), 
                                        .id = "traj", 
                                        deaths = NA, date = NA)) %>% 
                          return
                      }) %>% 
      rename(traj_set = X1) 
    
  } else{
    SEIR.sim <- simulate(covid_simulate,
                         nsim = nsim,
                         time = intervention.forecast@times,
                         format = "data.frame",
                         params = full_params)
  }
  
  
 SEIR.sim %<>% 
   mutate(paramset = variable_params[[i, "paramset"]],
          mif2_iter = variable_params[[i, "mif2_iter"]],
          date = round(as.Date(day, origin = variable_params[i, ]$sim_start)))
     
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



