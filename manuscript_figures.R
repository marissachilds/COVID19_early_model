# Set up ----
needed_packages <- c(
  "pomp"
  , "plyr"
  , "dplyr"
  , "ggplot2"
  , "magrittr"
  , "scales"
  , "lubridate"
  , "tidyr"
  , "gridExtra"
  , "data.table"
  , "grid"
  , "stringr")

lapply(needed_packages, library, character.only = TRUE)

# colors 
fig2_colors = c("#B40F20", "#E58601", "#046C9A", "grey30") # old colors
fig3_colors = c("#273046",
                "#366d8a",
                "#78B7C5")
fig4_colors = c("grey30",
                "#D67236",
                "#0b775e")

source("ggplot_theme.R")
# fitting        <- F
# fit.E0         <- F
# source("COVID_pomp.R")

# read deaths data used for fitting 
epi_df <- fread("https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-counties.csv") %>% 
  mutate(date = as.Date(date)) %>% 
  filter(county == "Santa Clara") %>% 
  filter(date <= as_date("2020-07-06")) %>%
  mutate(deaths_cum = deaths, cases_cum = cases) %>% 
  mutate(deaths = deaths_cum - lag(deaths_cum), cases = cases_cum - lag(cases_cum)) %>% 
  replace_na(list(deaths = 0))  %>%
  replace_na(list(cases = 0)) %>%
  mutate(raw_deaths = deaths, raw_cases = cases) #%>%

# Figure 3: R0 and Reff from different fit dates ----
all_fits_trajs = data.frame(
  fit_name = c("output/Santa_Clara_TRUE_FALSE_2020-04-01_2021-02-08_final.Rds",
               "output/Santa_Clara_TRUE_FALSE_2020-04-08_2021-02-05_final.Rds",
               "output/Santa_Clara_TRUE_FALSE_2020-04-15_2021-02-05_final.Rds",
               "output/Santa_Clara_TRUE_FALSE_2020-04-22_2021-02-05_final.Rds",
               "output/Santa_Clara_TRUE_FALSE_2020-04-29_2021-02-08_final.Rds",
               "output/Santa_Clara_TRUE_FALSE_2020-05-06_2021-02-05_final.Rds",
               "output/Santa_Clara_TRUE_FALSE_2020-05-13_2021-02-05_final.Rds",
               "output/Santa_Clara_TRUE_FALSE_2020-05-20_2021-02-05_final.Rds",
               "output/Santa_Clara_TRUE_FALSE_2020-05-27_2021-02-05_final.Rds",
               "output/Santa_Clara_TRUE_FALSE_2020-06-03_2021-02-05_final.Rds",
               "output/Santa_Clara_TRUE_FALSE_2020-06-10_2021-02-05_final.Rds",
               "output/Santa_Clara_TRUE_FALSE_2020-06-17_2021-02-05_final.Rds",
               "output/Santa_Clara_TRUE_FALSE_2020-06-24_2021-02-05_final.Rds")) %>% 
  separate(fit_name, into = c(NA, NA, NA, NA, "fit_date", NA), sep = "_", remove = FALSE) %>% 
  mutate(fit_date = as.Date(fit_date),
         traj_name = gsub(".Rds", "_filter_traj.Rds", fit_name))

variable_params_all <- all_fits_trajs %>% 
  mdply(function(fit_name, fit_date, traj_name){
    variable_params <- readRDS(fit_name)[["variable_params"]] 
  }) 

breaks_fun <- function(x) {
  if (max(x) > 20) { # breaks for cases
    seq(0, 300, by = 50)
  } else if(max(x) > 7) { # breaks for deaths
    seq(0, 10, by = 2)
  } else{ # breaks for reproduction number
    seq(0, 7, 1)
  }
}

figure1 <- variable_params_all %>% select(fit_date, R0, Reff, log_lik) %>%
  pivot_longer(R0:Reff) %>% 
  mutate(plot_panel = "Reproduction number\nestimates",
         date = fit_date) %>%
  group_by(fit_date) %>% 
  filter(log_lik > max(log_lik) - 2) %>% ungroup %>%
  {ggplot(data = .) +
  geom_violin(aes(x = fit_date, 
                  y = value, fill = name, 
                  color = name,
                  group = interaction(fit_date, name)),
              scale = "width", #"area", 
              # draw_quantiles = 0.5,
              # color = "black", 
              width = 5,
              trim = TRUE, 
              position = "identity") + 
  facet_grid(plot_panel ~ ., scales  = "free_y", 
             switch = "y", as.table = FALSE)  +
  geom_point(data = epi_df %>% 
               filter(date > as.Date("2020-03-01") & date < as.Date("2020-06-24")) %>% 
               select(date, cases, deaths) %>% pivot_longer(!date) %>% 
               mutate(plot_panel = stringr::str_to_title(name)),
             aes(x = date, y = value), color = "grey50") +
  geom_line(data = epi_df %>% 
              filter(date > as.Date("2020-03-01") & date < as.Date("2020-06-24")) %>% 
              select(date, cases, deaths) %>% pivot_longer(!date) %>% 
              mutate(plot_panel = stringr::str_to_title(name)),
            aes(x = date, y = value), color = "grey50") + 
  geom_hline(data = data.frame(plot_panel = "Reproduction number\nestimates", y = 1), 
             aes(yintercept = y), lty = "dashed") + 
      geom_vline(data = expand.grid(fit_date = pull(., fit_date) %>% unique,
                                   plot_panel = c("Deaths", "Cases")),
                 aes(xintercept = fit_date),
                 lty = "dotted", color = "grey60")+ 
  scale_color_manual(values = fig2_colors[1:2], name = "Estimate",
                     labels = c(expression(paste(R[0])),
                                expression(paste(R[eff]))),
                     aesthetics = c("colour", "fill")) + 
  scale_y_continuous(breaks = breaks_fun) +
  ylab("") +
  xlab("Last day of data used in fit") +
  theme(text=element_text(size=20),
        axis.title.y.right = element_text(hjust=0.8, vjust = 1),
        legend.text.align = 0,
        legend.position = c(0.08, 0.9),
        plot.margin = margin(t = 5.5, r = 5.5, b = 0, l = 5.5),
        strip.placement =  "outside", 
        strip.text = element_text(hjust = 0.5, vjust = 0.5)) }
  

# figure1

figure1 = ggplot_gtable(ggplot_build(figure1))
figure1$heights[7] = 2.2*figure1$heights[7]
# grid.draw(figure1)

ggsave("figures/R_ests.pdf", plot = figure1,
       device = "pdf", width = 15, height = 6,
       units = "in")

# SUPPLEMENT FIGURES ----

# dates to use for supplemental diagnostic plots
diag_dates <- as.Date(c("2020-04-01", "2020-05-13", "2020-06-24"))

# grab the mif traces 
mif_traces_all <- mdply(all_fits_trajs, function(fit_name, fit_date, traj_name){
  readRDS(fit_name)[["mif_traces"]] %>% 
    return 
}, .id = NULL) %>% 
  select(loglik, beta0, soc_dist_level_sip, E_init, 
         paramset, mif2_iter, mif_step, fit_date) 


# mif LL traces ----
{mif_traces_all  %>%
  filter(fit_date %in% diag_dates) %>%
  group_by(fit_date, paramset, mif2_iter) %>% 
  mutate(LL_max = max(loglik, na.rm = TRUE)) %>%
  group_by(fit_date) %>% 
  mutate(LL_max_std = LL_max - max(loglik, na.rm = TRUE)) %>%
  filter(paramset <= 10) %>%
  filter(mif_step > 30 & LL_max_std > -5) %>% 
  ggplot(aes(x = mif_step, 
             # y = asinh(loglik),
             y = loglik,
             group = interaction(fit_date, paramset, mif2_iter),
             color = as.factor(paramset))) +
  geom_line()  + 
  facet_wrap(~fit_date, scales = "free") +
  theme_bw()} %>% 
  ggsave(filename = "figures/mifs_LL_traces.pdf", width = 12, height = 6)

# violins of fit parameters ----
{variable_params_all %>% 
  select(-c(int_length2, int_start2, int_length1, iso_start)) %>% 
  group_by(fit_date) %>% 
  filter(log_lik > max(log_lik) -2 ) %>% 
  mutate_if(is.Date, ~as.numeric(difftime(.x, as.Date("2019-12-31")))) %>%
  # select(c(everything(), paramset, mif2_iter)) %>% colnames
  pivot_longer(cols = c(sim_start:beta0est, log_lik:soc_dist_level_sip, log_lik.se)) %>% 
  # filter(!(name %in% c("log_lik", "log_lik.se"))) %>%
  filter(name %in% c("beta0est", "E_init", "soc_dist_level_sip")) %>%
  mutate(name = ifelse(name == "beta0est", "beta0", name)) %>%
  # filter(name %in% c("E_init", "R0", "Reff")) %>%
  ggplot(aes(y = value, x = fit_date, 
             group = fit_date, fill = fit_date)) + 
  geom_violin() + 
  facet_wrap(~name, scales = "free") + 
  theme_bw()} %>% 
  ggsave("figures/fit_param_violins.pdf", plot =., width = 12, height = 4)

# mif replicates traces ----
{mif_traces_all  %>%
  group_by(fit_date, paramset, mif2_iter) %>% 
  mutate(LL_final = min(loglik, na.rm = TRUE)) %>%  
  ungroup %>% 
  pivot_longer(cols = loglik:E_init) %>% 
  filter(fit_date %in% diag_dates) %>%
  filter(name != "loglik") %>%
  filter(paramset <= 6) %>% 
  ggplot(aes(x = mif_step, y = value,
             color = as.factor(fit_date),
             group = interaction(fit_date, paramset, mif2_iter))) + 
  geom_line(alpha = 1) + 
  facet_grid(name ~ paramset, scales = "free_y") + 
  theme_bw()} %>% 
  ggsave("figures/mif_replicate_traces.pdf", plot = ., width = 12, height = 6)

# mif replicate LLs ??? ---- 
range_length = function(x){
  max(x) - min(x)
}
n_unique = function(x){
  length(unique(x))
}

variable_params_all %>% 
  filter(fit_date %in% diag_dates) %>% 
  select(-fit_name, - traj_name) %>% 
  group_by(fit_date, paramset) %>% 
  mutate_if(is.Date, ~as.numeric(difftime(.x, as.Date("2019-12-31")))) %>%
  select(-starts_with("int_length"), -iso_start, -int_start2, -log_lik.se,
         -c(beta0est, R0, Reff, E_init, soc_dist_level_sip)) %>% 
  select(fit_date, paramset, mif2_iter, log_lik, everything()) %>% 
  group_by(fit_date, paramset) %>% 
  summarise_at(vars(log_lik:mu), 
               list(range = range_length, mean = mean)) %>% 
  ungroup %>% 
  {select(., -c(summarise_all(., n_unique) %>%
                  unlist() %>% 
                  magrittr::equals(1) %>% 
                  which %>% 
                  names))} %>% 
  filter(fit_date == as.Date("2020-05-13")) %>% 
  select(-fit_date, -paramset) %>% 
  filter(log_lik_mean > max(log_lik_mean) - 10) %>% 
  pivot_longer(sim_start_mean:mu_mean) %>% 
  ggplot(aes(x = value, y = log_lik_range)) + 
  geom_point() + 
  facet_wrap(~name, scales = "free")

# profiles: parameter values vs LLs ----
{variable_params_all %>% 
  group_by(fit_date) %>% 
  mutate(LL_rel = log_lik - max(log_lik)) %>%
  filter(LL_rel > -2) %>%
  select(-c(int_length2, int_start2, int_length1, iso_start)) %>% 
  mutate_if(is.Date, ~as.numeric(difftime(.x, as.Date("2019-12-31")))) %>% 
  pivot_longer(cols = c(sim_start:beta0est, R0:soc_dist_level_sip)) %>% 
  filter(fit_date %in% diag_dates) %>%
  ggplot(aes(x = value, 
             # y = asinh(LL_rel), 
            y = LL_rel,
             color = as.factor(fit_date))) + 
  geom_point(alpha = 0.2) + 
  theme_bw() +
  facet_wrap(~name, scales = "free")} %>% 
  ggsave(filename = "figures/param_profiles.pdf", width = 12, height = 8)

# LL spread ----
{variable_params_all %>% 
  group_by(fit_date) %>% 
    mutate(LL_rel = log_lik - max(log_lik)) %>% 
    filter(LL_rel > -2) %>%
  ggplot(aes( x = -log_lik, y = log_lik.se,
              color = fit_date)) +
  geom_point() + 
  scale_y_continuous(trans = "log") + 
  scale_x_continuous(trans = "log") +
  theme_bw()} %>% 
  ggsave(filename = "figures/LL_spread_trunc.pdf", 
         plot = ., width = 10, height = 8)
  
# S over time ----
# warning: this takes a while to run 
pct_S_trajs <- data.frame(fit_name = c("output/Santa_Clara_TRUE_FALSE_2020-04-01_2021-02-08_final.Rds",
                                       "output/Santa_Clara_TRUE_FALSE_2020-05-06_2021-02-05_final.Rds", 
                                       "output/Santa_Clara_TRUE_FALSE_2020-06-10_2021-02-05_final.Rds"),
                          traj_name = c("output/Santa_Clara_TRUE_FALSE_2020-04-01_2021-02-08_final_filter_traj.Rds",
                                        "output/Santa_Clara_TRUE_FALSE_2020-05-06_2021-02-05_final_filter_traj.Rds",
                                        "output/Santa_Clara_TRUE_FALSE_2020-06-10_2021-02-05_final_filter_traj.Rds")) %>% 
  mdply(function(traj_name, fit_name){
    paramset_max <- 100000 # easy way to limit the number of paramsets to speed things up, 
    data_covar <- readRDS(fit_name)[["data_covar"]]
    variable_params <- readRDS(fit_name)[["variable_params"]]
    pop_size <- readRDS(fit_name)[["fixed_params"]] %>% 
      filter(param == "N") %>% pull(value) %>% first
    all_trajs <- readRDS(traj_name)
    fit_date <- as.Date(strsplit(fit_name, split = "_")[[1]][[5]])
    trajs <- variable_params %>% mutate(rowid = 1:nrow(.)) %>% 
      filter(log_lik > max(log_lik) - 2) %>% 
      arrange(desc(log_lik)) %>% 
      slice(1:paramset_max) %>% 
      select(rowid, log_lik) %>% 
      mdply(function(rowid, log_lik){
        all_trajs[[rowid]] %>% adply(3, function(traj){
          traj %>% t %>% data.frame(., day = as.numeric(rownames(.)))}) %>% 
          left_join(data_covar[[rowid]]$county_data %>% select(day, date)) %>%
          mutate(S_pct = S/pop_size) %>%
          group_by(date) %>% 
          summarise(mid = quantile(S_pct, 0.5),
                    lwr = quantile(S_pct, 0.025),
                    upr = quantile(S_pct, 0.975)) %>% 
          mutate(paramset = rowid,
                 log_lik  = log_lik,
                 fit_date = fit_date)}) %>% return
  })
 
{pct_S_trajs %>% 
  group_by(fit_date) %>% 
  mutate(LL_rel = log_lik - max(log_lik)) %>%
  ggplot(aes(x = date, y = mid, ymin = lwr, ymax = upr, 
             group = interaction(fit_date, paramset))) + 
  geom_ribbon(alpha = 0.2) +
  geom_violin(data = pct_S_trajs %>% 
                filter(date == fit_date),
              mapping = aes(date + as.difftime(18, units = "days"), 
                            y = mid, group = fit_date),
              adjust = 1) +
  geom_line(alpha = 0.1) +
  theme_bw() +
  # theme(legend.position = "none") +
  scale_fill_continuous(type = "viridis") +
  ylab("Percent S remaining") +
  facet_wrap(~fit_date, ncol = 1)} %>%
  {ggsave("figures/pct_S.pdf", ., width = 12, height = 7)}

