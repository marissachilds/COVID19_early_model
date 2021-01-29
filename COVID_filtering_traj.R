# Script to draw from the smoothing distribution from existing .Rds files with model fits
# Output is saved filtering trajectories

# required parameters, only set them if the file isn't being sourced
if(sys.nframe() == 0L){
  n.particles = 1000
  n.filter.traj = 10
  rds.name = "output/Santa Clara_TRUE_FALSE_2020-04-01_2021-01-26_final.Rds"
  seed.val = 100001
}

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

# Load the saved fit 
prev.fit           <- readRDS(rds.name)
variable_params    <- prev.fit[["variable_params"]] %>% 
  mutate(rowid = 1:n()) # add a row id to be able to match up with the data and covariate tables stored in data_covar
fixed_params       <- prev.fit[["fixed_params"]]
dynamics_summary   <- prev.fit[["dynamics_summary"]]
mif_traces         <- prev.fit[["mif_traces"]] 
data_covar         <- prev.fit[["data_covar"]]
pomp_components    <- prev.fit[["pomp_components"]]

set.seed(seed.val)
# loop over the rows of variable params
traj.all <- foreach(i = 1:nrow(variable_params), .combine = list, .multicombine = TRUE) %dopar%  {

    # combined the fixed params and variable params into a vector
  full_params = c(fixed_params %>% 
                    filter(paramset == variable_params[[i, "paramset"]] & 
                             mif2_iter == variable_params[[i, "mif2_iter"]]) %>% 
                    select(param, value) %>% tibble::deframe(), 
                  variable_params[i,] %>% rename(beta0 = beta0est) %>% 
                    select(beta0, Ca, alpha, mu, delta, 
                           E_init, soc_dist_level_sip) %>% unlist)
  
  # create the pomp object
  covid_fitting <- do.call(pomp, 
                           args = c(pomp_components, 
                                    list(data = data_covar[[variable_params[i,]$rowid]]$county_data, 
                                         covar = data_covar[[variable_params[i,]$rowid]]$covar_table)))
  # pfilter to get the trajectories
  trajs <- replicate(n.filter.traj, 
                     filter.traj(pfilter(covid_fitting, 
                                         params = full_params,
                                         filter.traj = TRUE, Np = n.particles))) 
  # returns array with dimensions number of states X 1 X time step (always includes 1 and then all of the covid_fitting@times) X n.filter.traj
  # for now, lets just combine them into a list, after dropping the 2nd dimension
  return(trajs[,1,,])
}

# save the filtering trajectories as the rds file name with _filter_traj appended to it
saveRDS(object = traj.all,
        file = strsplit(rds.name, ".", fixed = TRUE)[[1]] %>% 
          {paste0(c(.[-length(.)], "_filter_traj.", .[length(.)]), 
                  collapse = "")})


