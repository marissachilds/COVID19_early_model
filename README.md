# COVID19_early_model

Workflow 

1. Run fit.sh which uses the dates in fit_dates.csv and COVID_fit.R to fit the model weekly. Fits are saved in output/fits
2. Run filtering_traj.sh which uses fit_names.csv and COVID_filtering_traj.R to get filtering trajectories. Filtering trajectories are saved in output/filtering_trajectories.
3. Run assess.sh which uses sim-args.csv and sims-model_assess.R to run simulations for assessing model performace and saves them in output/simulations. 
4.  Run counterfactual_sims.R, inf_iso_sims.R, intervention_sims.R which each run simulations and save them in output/simulations.
5. Run manuscript_figures.R which creates figures from saved output and saves the figures in the figures directory.