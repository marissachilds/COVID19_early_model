# script for setting parameters for COVID_simulate 
####
## Parameters
####
seed.val           <- 10001
focal.county       <- "Santa Clara"
fitting            <- FALSE   ## Small change in pomp objects if fitting or simulating
## TRUE if COVID_fit previously run, FALSE if COVID_fit was just run and global environment is still full
use.rds            <- TRUE    
rds.name           <- "output/Santa Clara_TRUE_FALSE_2020-04-22_2021-01-08_final.Rds"
# more.params.uncer  <- FALSE   ## Fit with more (FALSE) or fewer (TRUE) point estimates for a number of parameters
# uncer.within.set   <- TRUE    ## Use parameter estimates from all mif2 runs for all parameter sets
nsim               <- 100       ## Number of epidemic simulations for each parameter set
fit.E0             <- TRUE      ## Was E0 also fit?

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
# sim_length       <- 400     ## How many days to run the simulation
state.plot         <- "D"     ## State variable for plotting (Hospit [H], Death [D], or Cases [C])

loglik.thresh      <- 2       ## Keep parameter sets with a likelihood within top X loglik units, to only fit with MLE, use 0
params.all         <- TRUE    ## Keep all fitted parameters above loglik thresh?...
nparams            <- 200     ## ...if FALSE, pick the top X by loglik to use
# plot.log10       <- FALSE   ## Plot on a log10 scale or not
