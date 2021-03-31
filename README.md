# COVID19_early_model

## Workflow 

1. Run fit.sh which uses the dates in fit_dates.csv and COVID_fit.R to fit the model weekly. Fits are saved in output/fits
2. Run filtering_traj.sh which uses fit_names.csv and COVID_filtering_traj.R to get filtering trajectories. Filtering trajectories are saved in output/filtering_trajectories.
3. Run assess.sh which uses sim-args.csv and sims-model_assess.R (which in turn uses COVID_simulate.R) to run simulations for assessing model performace and saves them in output/simulations. 
4.  Run counterfactual_sims.R, inf_iso_sims.R, intervention_sims.R (which in turn use COVID_simulate.R) which each run simulations and save them in output/simulations.
5. Run manuscript_figures.R which creates figures from saved output and saves the figures in the figures directory.

Saved output from steps 1 - 4 is available at https://drive.google.com/drive/folders/1mHu97K_FlPtcZVqlUPX5CKrrgKazT1N1?usp=sharing

## Other files
- COVID_pomp.R contains the model components as Csnippets.
- ggplot_theme.R sets a general ggplot theme. 
- location_params.csv contains Santa Clara County specific parameters that get used in COVID_fit.R
- params.csv defines the transmission parameters, including ranges for some parameters
- us-counties.txt is the data from the New York Times database: The New York Times. (2021). Coronavirus (Covid-19) Data in the United States. Retrieved from https://github.com/nytimes/covid-19-data.
- variable_params_less.R creates the Sobol sequences for sampled parameters and is used in COVID_fit.R