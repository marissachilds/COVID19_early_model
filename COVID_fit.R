##################################
## Fit COVID epidemic with pomp ##
##################################

args<-commandArgs(TRUE)

set.seed(879897)
fitting           <- TRUE          ## Small change in pomp objects if fitting or simulating
last_date         <- as.Date(args[1])  ## Last possible date to consider for this model
fit.E0            <- TRUE     ## Also fit initial # that starts the epidemic?
## more.params.uncer = FALSE is more supported, uses parameter ranges with more research and reacts to choice of focal.county if possible
## !!!!! For FALSE update parameters in location_params.csv
## more.params.uncer = TRUE  is less suppored, raw parameter values that can be adjusted manually
usable.cores      <- 2        ## Number of cores to use to fit
fit_to_sip        <- TRUE     ## Fit beta0 and shelter in place simultaneously?
import_cases      <- FALSE    ## Use importation of cases?
n.mif_runs        <- 6        ## mif2 fitting parameters
n.mif_length      <- 300
n.mif_particles   <- 1000
n.mif_rw.sd       <- 0.02
n.mif_particles_LL<- 10000     ## number of particles for calculating LL (10000 used in manuscript, 5000 suggested to debug/check code)
focal.county      <- "Santa Clara"  ## County to fit to
## !!! Curently parameters exist for Santa Clara, Miami-Dade, New York City, King, Los Angeles
## !!! But only Santa Clara explored
# county.N        <- 1.938e6         ## County population size
## !!! Now contained within location_params.csv
nparams           <- 200             ## number of parameter sobol samples (more = longer)
nsim              <- 200             ## number of simulations for each fitted beta0 for dynamics
download.new_data <- FALSE           ## Grab up-to-date data from NYT?

print(last_date)
print(class(last_date))

needed_packages <- c(
    "pomp"
  , "plyr"
  , "dplyr"
  , "ggplot2"
  , "magrittr"
  , "scales"
  , "lubridate"
  , "tidyr"
  , "foreach"
  , "doParallel"
  , "data.table")

lapply(needed_packages, require, character.only = TRUE)

## Be very careful here, adjust according to your machine
if(Sys.getenv('SLURM_JOB_ID') != ""){
  registerDoParallel(cores = (Sys.getenv("SLURM_NTASKS_PER_NODE")))
}else{
  registerDoParallel(cores = usable.cores)  
}

## Bring in pomp objects
source("COVID_pomp.R")
 
if (download.new_data) {
deaths <- fread("https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-counties.csv")    
} else {
deaths  <- read.csv("us-counties.txt")   
}
deaths    <- deaths %>% mutate(date = as.Date(date)) %>% filter(county == focal.county)
  
# last_date can't be later than the last date in existing data
last_date <- min(as.Date(last_date), max(deaths$date))
deaths    <- deaths %>% filter(date <= as.Date(last_date)) 
    
params <- read.csv("params.csv", stringsAsFactors = FALSE)
params <- params %>% mutate(Value = sapply(est, function(x) eval(parse(text = x))))

fixed_params        <- params$Value
names(fixed_params) <- params$Parameter
if (!import_cases) {fixed_params["import_rate"] <- 0}

location_params     <- read.csv("location_params.csv", stringsAsFactors = FALSE)
location_params     <- location_params %>% filter(location == focal.county)

fixed_params        <- c(fixed_params, N = location_params[location_params$Parameter == "N", ]$est)

## latest possible assumed epidemic start date
latest.sim_start <- as.Date(location_params[location_params$Parameter == "sim_start", ]$upr, origin = "2019-12-31")

source("variable_params_less.R")

## Repeat each entry for the number of repeat mif2
variable_params <- variable_params[rep(seq_len(nrow(variable_params)), each = n.mif_runs), ] %>% 
  mutate(mif2_iter = rep(seq(1, n.mif_runs), nparams))

## Run parameters
sim_start  <- variable_params$sim_start
sim_length <- 500
sim_end    <- sim_start + sim_length

## containers for results
SEIR.sim.ss.t.ci <- data.frame(
  name     = character(0)
, lwr      = numeric(0)
, est      = numeric(0)
, upr      = numeric(0)
, paramset = numeric(0))

fit.out <- foreach(i = 1:nrow(variable_params)) %dopar%  {
    
library(pomp)
library(dplyr)

  ## Adjust the data for the current start date
county.data <- deaths %>% 
  mutate(day = as.numeric(date - variable_params[i, ]$sim_start)) %>% 
  filter(day > 0) %>%
  dplyr::select(day, date, deaths) %>% 
  mutate(deaths_cum = deaths) %>% 
  mutate(deaths = deaths_cum - lag(deaths_cum)) %>% 
  replace_na(list(deaths = 0)) %>%
  dplyr::select(-deaths_cum) %>% 
  mutate(deaths = ifelse(date <= latest.sim_start, NA, deaths))

## Create intervention covariate table for the full forecast
if (fit_to_sip) {
  intervention.forecast <- with(variable_params[i, ], {

 covariate_table(
  day              = 1:sim_length
  
, intervention     = c(
      # No intervention until intervention start time
    rep(0, int_start1 - sim_start)                   
      # Intervention style 1
  , rep(1, int_length1)
      # Intervention style 2
  , rep(2, sim_length - (int_start2 - sim_start))
      # Post intervention close
)
  
, isolation        = rep(NA, sim_length)
, iso_severe_level = rep(NA, sim_length)      # % of contats that severe cases maintain
, iso_mild_level   = rep(NA, sim_length)   # % of contats that mild cases maintain

, soc_dist_level_wfh = rep(soc_dist_level_wfh, sim_length) 

, thresh_H_start   = rep(NA, sim_length)
, thresh_H_end     = rep(NA, sim_length)

, thresh_int_level = rep(NA, sim_length)
, back_int_level   = rep(NA, sim_length)
  
, order            = "constant"
, times            = "day"
  )

})
} else {
  intervention.forecast <- with(variable_params[i, ], {

 covariate_table(
  day              = 1:sim_length
  
, intervention     = c(
      # No intervention until intervention start time
    rep(0, int_start1 - sim_start)                   
      # Intervention style 1
  , rep(1, int_length1)
      # Intervention style 2
  , rep(2, sim_length - (int_start2 - sim_start))
      # Post intervention close
)
  
, isolation        = rep(NA, sim_length)
, iso_severe_level = rep(NA, sim_length)      # % of contats that severe cases maintain
, iso_mild_level   = rep(NA, sim_length)   # % of contats that mild cases maintain

, soc_dist_level_wfh = rep(soc_dist_level_wfh, sim_length) 
, soc_dist_level_sip = rep(soc_dist_level_sip, sim_length)

, thresh_H_start   = rep(NA, sim_length)
, thresh_H_end     = rep(NA, sim_length)

, thresh_int_level = rep(NA, sim_length)
, back_int_level   = rep(NA, sim_length)
  
, order            = "constant"
, times            = "day"
  )

})  
}

covid.fitting <- county.data %>%
  pomp(
    time       = "day"
  , t0         = 1
  , covar      = intervention.forecast
  , rprocess   = euler(sir_step, delta.t = 1/6)
  , rmeasure   = rmeas_deaths
  , dmeasure   = dmeas_deaths 
  , rinit      = sir_init
  , partrans   = par_trans
  , accumvars  = accum_names
  , paramnames = param_names
  , statenames = state_names) 

if (variable_params[i, ]$beta0est == 0) {

if (!more.params.uncer) {

mifs_temp <- try(silent = TRUE, { covid.fitting %>%
  mif2(
    t0      = 1
  , params  = c(
    c(fixed_params
    , Ca    = variable_params[i, ]$Ca
    , alpha = variable_params[i, ]$alpha
    , delta = variable_params[i, ]$delta
    , mu    = variable_params[i, ]$mu
      )
  , {
    if (fit_to_sip) {
      if (fit.E0) {
    c(beta0              = rlnorm(1, log(0.7), 0.17)
    , soc_dist_level_sip = rlnorm(1, log(0.2), 0.2)
    , E_init             = runif(1, 0, 5))
      } else {
    c(beta0              = rlnorm(1, log(0.7), 0.17)
    , soc_dist_level_sip = rlnorm(1, log(0.2), 0.2)
    , E0                 = variable_params[i, ]$E0)        
      }
    } else {
      if (fit.E0) {
    c(beta0              = rlnorm(1, log(0.7), 0.17)
      , E_init             = runif(1, 0, 5))
      } else {
    c(beta0              = rlnorm(1, log(0.7), 0.17)
    , E0                 = variable_params[i, ]$E0)        
      }
    }
  }
  )
  , Np                  = n.mif_particles
  , Nmif                = n.mif_length
  , cooling.fraction.50 = 0.5
  , rw.sd               = {
    if (fit_to_sip) {
      if (fit.E0) {
    rw.sd(
      beta0              = n.mif_rw.sd
    , soc_dist_level_sip = n.mif_rw.sd
    , E_init             = ivp(n.mif_rw.sd)
      )       
      } else {
    rw.sd(beta0 = n.mif_rw.sd, soc_dist_level_sip = n.mif_rw.sd)        
      }
    } else {
      if (fit.E0) {
    rw.sd(
      beta0  = n.mif_rw.sd
    , E_init = ivp(n.mif_rw.sd)
      )        
      } else {
    rw.sd(beta0 = n.mif_rw.sd)        
      }
    }
  }
)
  })
  
} else {
  
mifs_temp <- try(silent = TRUE, { covid.fitting %>% mif2(
    t0     = 1
  , params = c(
    c(fixed_params
      , Ca       = variable_params[i, ]$Ca
      , alpha    = variable_params[i, ]$alpha
      , delta    = variable_params[i, ]$delta
      , mu       = variable_params[i, ]$mu
      , lambda_a = variable_params[i, ]$lambda_a
      , lambda_s = variable_params[i, ]$lambda_s
      , lambda_m = variable_params[i, ]$lambda_m  
      )
  , {
    if (fit_to_sip) {
      if (fit.E0) {
    c(beta0              = rlnorm(1, log(0.7), 0.17)
    , soc_dist_level_sip = rlnorm(1, log(0.2), 0.2)
    , E_init             = runif(1, 0, 5))
      } else {
    c(beta0              = rlnorm(1, log(0.7), 0.17)
    , soc_dist_level_sip = rlnorm(1, log(0.2), 0.2)
    , E0                 = variable_params[i, ]$E0)        
      }
    } else {
      if (fit.E0) {
    c(beta0              = rlnorm(1, log(0.7), 0.17)
      , E_init             = runif(1, 0, 5))
      } else {
    c(beta0              = rlnorm(1, log(0.7), 0.17)
    , E0                 = variable_params[i, ]$E0)        
      }
    }
  }
  )
  , Np                  = n.mif_particles
  , Nmif                = n.mif_length
  , cooling.fraction.50 = 0.5
  , rw.sd               = {
    if (fit_to_sip) {
      if (fit.E0) {
    rw.sd(
      beta0              = n.mif_rw.sd
    , soc_dist_level_sip = n.mif_rw.sd
    , E_init             = ivp(n.mif_rw.sd))        
      } else {
    rw.sd(beta0 = n.mif_rw.sd, soc_dist_level_sip = n.mif_rw.sd)        
      }
    } else {
      if (fit.E0) {
    rw.sd(
      beta0  = n.mif_rw.sd
    , E_init = ivp(n.mif_rw.sd)
      )        
      } else {
    rw.sd(beta0 = n.mif_rw.sd)        
      }
    }
  }
)})

}
  
}

testforerror <- sum(class(mifs_temp) == "try-error")

if (testforerror == 0) {
  
mifs.ll <- try(silent = TRUE, {replicate(10, mifs_temp %>% pfilter(Np = n.mif_particles_LL) %>% logLik())})
mifs.ll <- try(silent = TRUE, {logmeanexp(mifs.ll, se = TRUE)})

mifs_temp.coef <- coef(mifs_temp)

if (fit_to_sip) {
 variable_params[i, "soc_dist_level_sip"] <- mifs_temp.coef["soc_dist_level_sip"]
}
if (fit.E0) {
variable_params[i, "E_init"]              <- mifs_temp.coef["E_init"]
}

variable_params[i, "beta0est"]   <- mifs_temp.coef["beta0"]
variable_params[i, "log_lik"]    <- mifs.ll[1]
variable_params[i, "log_lik.se"] <- mifs.ll[2]

mif_traces    <- traces(mifs_temp)

SEIR.sim <- do.call(
  pomp::simulate
  , list(
    object         = covid.fitting
    , times        = intervention.forecast@times
    , params = {
    if (!more.params.uncer) {
      c(fixed_params, c(
      beta0              = variable_params[i, "beta0est"]
    , soc_dist_level_sip = variable_params[i, "soc_dist_level_sip"]
    , Ca                 = variable_params[i, ]$Ca
    , alpha              = variable_params[i, ]$alpha
    , delta              = variable_params[i, ]$delta
    , mu                 = variable_params[i, ]$mu
      )
      , if (fit.E0) { 
          c(E_init = variable_params[i, ]$E_init)
        } else {
          c(E0     = variable_params[i, ]$E0)     
        }
        )
      } else {
      c(fixed_params, c(
      beta0              = variable_params[i, "beta0est"]
    , soc_dist_level_sip = variable_params[i, "soc_dist_level_sip"]
      , Ca               = variable_params[i, ]$Ca
      , alpha            = variable_params[i, ]$alpha
      , lambda_a         = variable_params[i, ]$lambda_a
      , lambda_s         = variable_params[i, ]$lambda_s
      , lambda_m         = variable_params[i, ]$lambda_m 
      )
      , if (fit.E0) { 
          c(E_init = variable_params[i, ]$E_init)
        } else {
          c(E0     = variable_params[i, ]$E0)     
        }
        ) 
      }}
    , nsim         = nsim
    , format       = "d"
    , include.data = F
    , seed         = 1001)) %>% {
      rbind(.,
         group_by(., day) %>%
           dplyr::select(-.id) %>%
           summarise_all(median) %>%
                    mutate(.id = "median"))
    }


## summarize epidemic
{
SEIR.sim.s  <- SEIR.sim  %>% 
  dplyr::group_by(.id) %>% 
  dplyr::summarize(total_H = sum(H))

SEIR.sim    <- left_join(SEIR.sim, SEIR.sim.s, by = ".id")

SEIR.sim    <- SEIR.sim %>% 
  dplyr::filter(
    total_H > 10
  ) %>% droplevels()

## Calc the summary statistics
SEIR.sim.ss <- SEIR.sim %>% 
  dplyr::filter(.id != "median") %>%
  dplyr::mutate(week   = day %/% 7) %>%
  dplyr::group_by(week, .id) %>%
  dplyr::summarize(
    ## Mean in hospitalizations by week
    sum_H = sum(H_new)
    ) %>%
  dplyr::group_by(.id) %>%  
  dplyr::mutate(
    ## Difference in hospitalizations at the level of the week
    diff_H = c(0, diff(sum_H))
    ) %>%
  dplyr::group_by(.id)

SEIR.sim.ss.1 <- SEIR.sim.ss %>% 
  dplyr::summarize(
    ## Maximum hospitalizations reached, summarized at the level of the week
    when_max_H = week[which.max(sum_H)]
    ## How many hospitalizations are reached in that week
  , max_H      = max(sum_H)
    ## First week we see a reduction in the number of hospitalizations from a runs _global_ rate peak
  , when_red_H = week[min(which(diff_H == max(diff_H)))]
    )

## Number S to calc R_eff
SEIR.sim.ss.2 <- SEIR.sim %>% 
  dplyr::filter(day == max(county.data$day)) %>%
  dplyr::select(.id, S, D)

names(SEIR.sim.ss.2)[c(2, 3)] <- c("S_now", "D_now")

SEIR.sim.ss.f <- left_join(SEIR.sim.ss.1, SEIR.sim.ss.2, by = ".id")

SEIR.sim.ss.3 <- SEIR.sim %>%
  filter(.id != "median") %>%
  dplyr::group_by(.id) %>% 
  dplyr::summarize(
    total_D = max(D)
  , total_R = max(R)
      )

SEIR.sim.ss.t   <- left_join(SEIR.sim.ss.f, SEIR.sim.ss.3, by = ".id")

## Get their CI
SEIR.sim.ss.t.s <- SEIR.sim.ss.t %>%
  pivot_longer(cols = -.id) %>%
  dplyr::group_by(name) %>%
  dplyr::summarize(
    lwr = quantile(value, c(0.025))
  , est = quantile(value, c(0.500))
  , upr = quantile(value, c(0.975))
  ) %>% 
  dplyr::mutate(paramset = variable_params[i, ]$paramset)

SEIR.sim.ss.t.s <- rbind(SEIR.sim.ss.t.s
  , data.frame(
      name = "prop_ext"
    , lwr  = nsim - nrow(SEIR.sim.ss.t)
    , est  = nsim - nrow(SEIR.sim.ss.t)
    , upr  = nsim - nrow(SEIR.sim.ss.t)
    , paramset = variable_params[i, ]$paramset))

SEIR.sim.ss.t.ci <- rbind(SEIR.sim.ss.t.ci, SEIR.sim.ss.t.s)
}

## Update variable params with R0 estimate
variable_params[i, ]$Reff <- with(variable_params[i, ], covid_R0(
  beta0est = beta0est, fixed_params = c(fixed_params, unlist(variable_params[i, ]))
  , sd_strength = soc_dist_level_sip
, prop_S = unlist(SEIR.sim.ss.t.s[SEIR.sim.ss.t.s$name == "S_now", 3]) / 
    (location_params[location_params$Parameter == "N", ]$est - 
        SEIR.sim.ss.t.s[SEIR.sim.ss.t.s$name == "total_D", 3])))

variable_params[i, ]$R0 <- with(variable_params[i, ], covid_R0(
  beta0est = beta0est, fixed_params = c(fixed_params, unlist(variable_params[i, ]))
  , sd_strength = 1, prop_S = 1))

}

## Clean up the objects to return 
fixed_params <- melt(fixed_params) %>% mutate(
  param      = names(fixed_params)
  , paramset  = variable_params$paramset[i]
  , mif2_iter = variable_params$mif2_iter[i]
)
SEIR.sim.ss.t.ci <- SEIR.sim.ss.t.ci %>%
  mutate(
    mif2_iter = variable_params$mif2_iter[i] 
  )
mif_traces <- as.data.frame(mif_traces) %>% mutate(
  paramset  = variable_params$paramset[i]
  , mif2_iter = variable_params$mif2_iter[i]
  , mif_step  = seq(1, nrow(mif_traces))
)

return(
list(
    variable_params  = variable_params[i, ]
  , fixed_params     = fixed_params
  , dynamics_summary = SEIR.sim.ss.t.ci
  , mif_traces       = mif_traces
  , data_covar       = list(county_data      = county.data
                            , covar_table      = intervention.forecast)
   )
)

}

variable_params.all  <- sapply(fit.out, FUN = function(x) x[1]) %>% rbindlist()
fixed_params.all     <- sapply(fit.out, FUN = function(x) x[2]) %>% rbindlist()
dynamics_summary.all <- sapply(fit.out, FUN = function(x) x[3]) %>% rbindlist()
mif_traces.all       <- sapply(fit.out, FUN = function(x) x[4]) %>% rbindlist()
data_covar.all       <- sapply(fit.out, FUN = function(x) x$data_covar, simplify = FALSE)

saveRDS(
   list(
    variable_params  = variable_params.all
  , fixed_params     = fixed_params.all
  , dynamics_summary = dynamics_summary.all
  , mif_traces       = mif_traces.all
  , data_covar       = data_covar.all
  , pomp_components  = list(time       = "day"
                            , t0         = 1
                            , rprocess   = euler(sir_step, delta.t = 1/6)
                            , rmeasure   = rmeas_deaths
                            , dmeasure   = dmeas_deaths 
                            , rinit      = sir_init
                            , partrans   = par_trans
                            , accumvars  = accum_names
                            , paramnames = param_names
                            , statenames = state_names)
   ), paste(
     paste("output/fits/"
           , paste(gsub(" ", "_", focal.county), fit_to_sip, more.params.uncer, last_date, Sys.Date(), "final", sep = "_")
         , sep = "")
     , "Rds", sep = "."))
 
