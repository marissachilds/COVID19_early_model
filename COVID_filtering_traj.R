# Script to draw from the smoothing distribution from existing .Rds files with model fits
# Output is saved filtering trajectories
needed_packages <- c(
  "pomp"
  , "plyr"
  , "dplyr"
  , "magrittr"
  , "tidyr"
  , "foreach"
  , "data.table"
  , "doParallel")

lapply(needed_packages, require, character.only = TRUE)

args<-commandArgs(TRUE)

n.particles <- 5000
n.filter.traj <- 100
rds.name <- paste0("output/", args[1])
usable.cores <- 1
print(rds.name)

if(Sys.getenv('SLURM_JOB_ID') != ""){
  registerDoParallel(cores = (Sys.getenv("SLURM_NTASKS_PER_NODE")))
}else{
  registerDoParallel(cores = usable.cores)  
}

# Load the saved fit 
prev.fit           <- readRDS(rds.name)
variable_params    <- prev.fit[["variable_params"]] %>% 
  mutate(rowid = 1:n()) # add a row id to be able to match up with the data and covariate tables stored in data_covar
fixed_params       <- prev.fit[["fixed_params"]]
dynamics_summary   <- prev.fit[["dynamics_summary"]]
mif_traces         <- prev.fit[["mif_traces"]] 
data_covar         <- prev.fit[["data_covar"]]
pomp_components    <- prev.fit[["pomp_components"]]

set.seed(100001)

# loop over the rows of variable params
traj.all <- foreach(
  i = 1:nrow(variable_params)) %dopar%  {
    tic <- Sys.time()                  
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
  toc <- Sys.time() 
  print(toc - tic)
  # returns array with dimensions number of states X 1 X time step (always includes 1 and then all of the covid_fitting@times) X n.filter.traj
  # for now, lets just combine them into a list, after dropping the 2nd dimension
  return(trajs[,1,,])
}



# save the filtering trajectories as the rds file name with _filter_traj appended to it
saveRDS(object = traj.all,
        file = strsplit(rds.name, ".", fixed = TRUE)[[1]] %>% 
          {paste0(c(.[-length(.)], "_filter_traj.", .[length(.)]), 
                  collapse = "")})


