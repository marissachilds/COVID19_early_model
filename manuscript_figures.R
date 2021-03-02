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

# read deaths data used for fitting 
epi_df <- fread("https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-counties.csv") %>% 
  mutate(date = as.Date(date)) %>% 
  filter(county == "Santa Clara") %>% 
  # filter(date <= as_date("2020-07-06")) %>%
  mutate(deaths_cum = deaths, cases_cum = cases) %>% 
  mutate(deaths = deaths_cum - lag(deaths_cum), cases = cases_cum - lag(cases_cum)) %>% 
  replace_na(list(deaths = 0))  %>%
  replace_na(list(cases = 0)) %>%
  mutate(raw_deaths = deaths, raw_cases = cases) #%>%

# figure 3: reproduction numbers ----
all_fits_trajs = data.frame(
  fit_name = c("output/fits/Santa_Clara_TRUE_FALSE_2020-04-01_2021-02-08_final.Rds",
               "output/fits/Santa_Clara_TRUE_FALSE_2020-04-08_2021-02-05_final.Rds",
               "output/fits/Santa_Clara_TRUE_FALSE_2020-04-15_2021-02-05_final.Rds",
               "output/fits/Santa_Clara_TRUE_FALSE_2020-04-22_2021-02-05_final.Rds",
               "output/fits/Santa_Clara_TRUE_FALSE_2020-04-29_2021-02-08_final.Rds",
               "output/fits/Santa_Clara_TRUE_FALSE_2020-05-06_2021-02-05_final.Rds",
               "output/fits/Santa_Clara_TRUE_FALSE_2020-05-13_2021-02-05_final.Rds",
               "output/fits/Santa_Clara_TRUE_FALSE_2020-05-20_2021-02-05_final.Rds",
               "output/fits/Santa_Clara_TRUE_FALSE_2020-05-27_2021-02-05_final.Rds",
               "output/fits/Santa_Clara_TRUE_FALSE_2020-06-03_2021-02-05_final.Rds",
               "output/fits/Santa_Clara_TRUE_FALSE_2020-06-10_2021-02-05_final.Rds",
               "output/fits/Santa_Clara_TRUE_FALSE_2020-06-17_2021-02-05_final.Rds",
               "output/fits/Santa_Clara_TRUE_FALSE_2020-06-24_2021-02-05_final.Rds")) %>% 
  separate(fit_name, into = c(NA, NA, NA, NA, "fit_date", NA), sep = "_", remove = FALSE) %>% 
  mutate(fit_date = as.Date(fit_date),
         traj_name = gsub("fits", "filter_trajectories", gsub(".Rds", "_filter_traj.Rds", fit_name)))

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

r_ests <- variable_params_all %>% select(fit_date, R0, Reff, log_lik) %>%
  pivot_longer(R0:Reff) %>% 
  mutate(plot_panel = "Reproduction number\nestimates",
         date = fit_date) %>%
  group_by(fit_date) %>% 
  filter(log_lik > max(log_lik) - 2) %>% ungroup %>%
  {ggplot(data = .) +
      geom_violin(aes(x = fit_date, 
                      y = value, fill = name, 
                      group = interaction(fit_date, name)),
                  scale = "width", #"area", 
                  color = "black",
                  width = 5,
                  size = 0.5,
                  trim = TRUE, 
                  position = "identity") + 
      facet_grid(plot_panel ~ ., scales  = "free_y", 
                 switch = "y", as.table = FALSE)  +
      geom_point(data = epi_df %>% 
                   filter(date > as.Date("2020-03-01") & date <= as.Date("2020-06-24")) %>% 
                   select(date, cases, deaths) %>% pivot_longer(!date) %>% 
                   mutate(plot_panel = ifelse(name == "deaths", "Daily deaths", "Daily reported\ncases")),
                 aes(x = date, y = value), color = "grey50") +
      geom_line(data = epi_df %>% 
                  filter(date > as.Date("2020-03-01") & date <= as.Date("2020-06-24")) %>% 
                  select(date, cases, deaths) %>% pivot_longer(!date) %>% 
                  mutate(plot_panel = ifelse(name == "deaths", "Daily deaths", "Daily reported\ncases")),
                aes(x = date, y = value), color = "grey50") + 
      geom_hline(data = data.frame(plot_panel = "Reproduction number\nestimates", y = 1), 
                 aes(yintercept = y), lty = "dashed") + 
      geom_vline(data = expand.grid(fit_date = pull(., fit_date) %>% unique,
                                    plot_panel = c("Daily deaths", "Daily reported\ncases")),
                 aes(xintercept = fit_date),
                 lty = "dotted", color = "grey60")+ 
      scale_color_manual(values = alpha(fig2_colors[1:2], 0.5), 
                         name = "Estimate",
                         labels = c(expression(paste(R[0])),
                                    expression(paste(R[eff]))),
                         aesthetics = c("colour", "fill")) + 
      scale_y_continuous(breaks = breaks_fun) +
      ylab("") +
      xlab("Last day of data used in fit") +
      theme(text=element_text(size=16),
            axis.title.y.right = element_text(hjust=0.8, vjust = 1),
            legend.text.align = 0,
            legend.position = c(0.08, 0.9),
            plot.margin = margin(t = 5.5, r = 5.5, b = 0, l = 5.5),
            strip.placement =  "outside", 
            strip.text = element_text(hjust = 0.5, vjust = 0.5)) }


r_ests

r_ests = ggplot_gtable(ggplot_build(r_ests))
r_ests$heights[7] = 2.2*r_ests$heights[7]
# grid.draw(r_ests)

ggsave("figures/R_ests.pdf", plot = r_ests,
       device = "pdf", width = 15, height = 6,
       units = "in")

# figure 4: model assessment for deaths ----
# due to memory limits, might need to read it in, only saving CIs/medians and error scores
assess_sims <- alply(list.files("./output/simulations", pattern = "2021-02-22" , full.names = T),
                     1,
                     function(sim_name){
                       print(sim_name)
                       summary_vars <- c("D_new", "I_new_sympt", "deaths")
                       sims <- readRDS(sim_name) %>% select(fit_date, paramset, mif2_iter, traj_set, .id, day, date, all_of(summary_vars))
                       # ci and median over all parameter sets
                       sims_ci <- sims %>%
                         group_by(fit_date, date) %>%
                         summarise_at(all_of(summary_vars), ~list(quantile(.x, c(0.025, 0.5, 0.975), na.rm = TRUE))) %>%
                         unnest(all_of(summary_vars)) %>%
                         mutate(q_names = paste0("q", gsub("%", "", names(D_new), fixed = T))) %>%
                         pivot_longer(all_of(summary_vars), names_to = "var")  %>%
                         pivot_wider(names_from = q_names) %>%
                         as.data.frame()
                       # medians for individual parameter sets
                       sims_mid <- sims %>%
                         group_by(fit_date, date, paramset, mif2_iter) %>%
                         summarise_at(all_of(summary_vars), median) %>%
                         as.data.frame()
                       # model eval
                       model_score <- sims %>%
                         mutate(days_since_fit = difftime(date, fit_date, units = "days")) %>%
                         filter(days_since_fit > 0 & days_since_fit <= 14) %>%
                         left_join(epi_df %>% select(date, deaths) %>% rename(obs_deaths = deaths),
                                   by = "date") %>%
                         rowwise() %>%
                         mutate(q_score = -2*dpois(obs_deaths, deaths) +
                             besselI(x = 2*deaths, nu = 0, expon.scaled = TRUE),
                           abs_error = abs(obs_deaths - deaths)) %>%
                         # calculate average from each simulation or each parameter set? currentl each simulation
                         group_by(fit_date, paramset, mif2_iter, traj_set, .id) %>% # 
                         summarise(q_score_mean = mean(q_score),
                                   mae = mean(abs_error),
                                   pct_over = sum(deaths > obs_deaths)/length(deaths),
                                   .groups = "drop")
                       return(list(ci = sims_ci,
                                   mid = sims_mid,
                                   score = model_score))
                     })

# # separate into different components 
assess_sims <- list(ci = ldply(assess_sims, function(l){l[["ci"]]}),
                  mid = ldply(assess_sims, function(l){l[["mid"]]}),
                  score = ldply(assess_sims, function(l){l[["score"]]}))

# saveRDS(assess_sims, "output/asses_sims_summarised_by_sim.Rds")
# readRDS("output/assess_sims_summarised_by_sim.Rds")

traj_dates <- as.Date(c("2020-04-01", "2020-04-22", "2020-05-13", 
                        "2020-06-03", "2020-06-24"))

deaths_traj_fig <- assess_sims[["ci"]] %>% 
  filter(fit_date %in% traj_dates) %>% 
  rowwise %>% 
  mutate(q97.5 = min(q97.5, 15)) %>%
  # filter(date <= max(fit_date) + as.difftime(14, units = "days")) %>%
  {ggplot(data = .) + 
      facet_wrap(~fit_date, nrow = 1, labeller = function(x){format(x, format = "%B %e")}) +
      # central 95% interval among all future sims
      geom_ribbon(mapping = aes(x = date, ymin = q2.5, ymax = q97.5),
                  data = filter(., date > fit_date & var == "deaths"),
                  alpha = 0.5, fill = "grey", color = NA) +
      # central 95% interval among all past trajectories
      geom_ribbon(mapping = aes(x = date, ymin = q2.5, ymax = q97.5),
                  data = filter(., date <= fit_date & var == "D_new"),
                  alpha = 0.5, fill = "grey", color = NA) +
      # median lines for individual paramsets among future sims
      geom_line(mapping = aes(x = date, y = deaths,
                              group = interaction(fit_date, paramset, mif2_iter)),
                data = assess_sims[["mid"]] %>%
                  filter(date > fit_date) %>%
                  filter(fit_date %in% traj_dates),
                alpha = 0.1, color = "grey30") +
      # median lines for individual paramsets among past trajectories
      geom_line(mapping = aes(x = date, y = D_new,
                              group = interaction(fit_date, paramset, mif2_iter)),
                data = assess_sims[["mid"]] %>%
                  filter(date <= fit_date) %>%
                  filter(fit_date %in% traj_dates),
                alpha = 0.01, color = "grey30")  +
      # add fit date lines
      geom_vline(data = select(., fit_date) %>% unique,
                 mapping = aes(xintercept = fit_date),
                 linetype = "dashed") + 
      # add data points
      geom_point(data.frame(fit_date = unique(.$fit_date)) %>% 
                   left_join(epi_df,
                             by = character()) %>% 
                   filter(date <= as.Date("2020-07-31")) %>%
                   mutate(diff_date = difftime(date, fit_date, units = "days"),
                          type = case_when(diff_date <= 0 ~ "pre",
                                           diff_date > 14 ~ "post", 
                                           T ~ "assess"),
                          type = factor(type, levels = c("pre", "assess", "post"))),
                 mapping = aes(x = date, y = deaths, color = type)) + 
      scale_colour_manual(values = c("black",# "#222222", #"#00A08A", # deaths past for fitting
                                     "#CA0020", # future predictions within 2 weeks of fit
                                     alpha("blue", 0.3)), # "#F4A582" #  #"#EF8A62", #"#F4A582", #alpha("#FF0000", 0.4),  # future predictions > 2 weeks of fit
                          labels = c("Deaths used for\nmodel fitting", "Future predicted,\ncompared",
                                     "Future predicted,\nuncompared"),
                          name = "") + 
      xlab("") + ylab("Daily deaths") + 
      ylim(0, 15) +
      scale_x_date(breaks = "2 months",
                   limits = as.Date(c("2020-02-15", NA)),
                   date_labels = "%b")}

# deaths_traj_fig

score_fig <- assess_sims[["score"]] %>%
  group_by(fit_date) %>% 
  mutate(pct_over_mean = mean(pct_over)) %>%
  {ggplot(data = .) + 
      geom_rect(data = select(., fit_date) %>% unique %>% 
                  filter(fit_date %in% traj_dates) %>%
                  mutate(date_min = fit_date - as.difftime(3.5, units = "days"), 
                         date_max = fit_date + as.difftime(3.5, units = "days")),
                mapping = aes(xmin = date_min, xmax = date_max),
                ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.5) +
      geom_violin(mapping = aes(x = fit_date, y = q_score_mean, group = fit_date,
                                size = as.numeric(fit_date %in% traj_dates),
                                fill = pct_over_mean),
                  color = "black", width = 7, scale = "area", adjust = 2) + 
      geom_text(data = summarise(., min_score = min(q_score_mean)) %>% 
                  filter(fit_date %in% traj_dates), 
                mapping = aes(x = fit_date, y = min_score - 0.1, 
                              label = format(fit_date, "%b %e")),
                inherit.aes = FALSE) +
      scale_size_binned(range = c(0.01, 0.8), guide = FALSE) +
      scale_fill_gradient2(low = muted("blue"), mid = "white", high = muted("red"),
                           midpoint = 0.5, guide = "colourbar", 
                           name = "Percent\nover predicted",
                           # breaks = c(0, 0.25, 0.5, 0.75),
                           labels = label_percent(accuracy = 1)) +
      xlab("Fitting Date") + ylab("Quadratic score")}

# score_fig

ggpubr::ggarrange(deaths_traj_fig, score_fig, nrow = 2, 
                  heights = c(6, 4)) %>% 
  ggsave(filename = "figures/deaths_assess.pdf", width = 16, height = 7)
# https://stats.stackexchange.com/questions/162429/sum-of-squared-poisson-probability-masses
# https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.67.3696&rep=rep1&type=pdf
# p_2_dict <- data.frame(lambda = 0:100) %>% rowwise() %>% 
#   mutate(approx = sum(dpois(0:1000, lambda)^2),
#          analytic = besselI(x = 2*lambda, nu = 0, expon.scaled = TRUE),
#          )

# figure 5: model assessment for cases ----

sims_lag <- 7
cases_scaling <- 25
cases_traj_fig <- assess_sims[["ci"]] %>% 
  filter(fit_date %in% traj_dates) %>% 
  mutate(date = date + as.difftime(sims_lag, units = "days")) %>%
  {ggplot(data = .) + 
      facet_wrap(~fit_date, nrow = 1, labeller = function(x){format(x, format = "%B %e")}) +
      # central 95% interval among all future sims
      geom_ribbon(mapping = aes(x = date, ymin = q2.5, ymax = q97.5),
                  data = filter(., date > fit_date & var == "I_new_sympt"),
                  alpha = 0.5, fill = "grey", color = NA) +
      # central 95% interval among all past trajectories
      geom_ribbon(mapping = aes(x = date, ymin = q2.5, ymax = q97.5),
                  data = filter(., date <= fit_date & var == "I_new_sympt"),
                  alpha = 0.5, fill = "grey", color = NA) +
      # median lines for individual paramsets among future sims
      geom_line(mapping = aes(x = date, y = I_new_sympt, 
                              group = interaction(fit_date, paramset, mif2_iter)),
                data = assess_sims[["mid"]] %>% 
                  mutate(date = date + as.difftime(sims_lag, units = "days")) %>%
                  filter(date > fit_date) %>%
                  filter(fit_date %in% traj_dates),
                alpha = 0.05) + 
      # median lines for individual paramsets among past trajectories
      geom_line(mapping = aes(x = date, y = I_new_sympt, 
                              group = interaction(fit_date, paramset, mif2_iter)),
                data = assess_sims[["mid"]] %>% 
                  mutate(date = date + as.difftime(sims_lag, units = "days")) %>%
                  filter(date <= fit_date) %>%
                  filter(fit_date %in% traj_dates),
                alpha = 0.05)  + 
      # add fit date lines
      geom_vline(data = select(., fit_date) %>% unique,
                 mapping = aes(xintercept = fit_date),
                 linetype = "dashed") + 
      # add data points
      geom_point(data.frame(fit_date = unique(.$fit_date)) %>% 
                   left_join(epi_df,
                             by = character()) %>% 
                   filter(date <= as.Date("2020-07-31")) %>%
                   mutate(diff_date = difftime(date, fit_date, units = "days"),
                          type = case_when(diff_date <= 0 ~ "pre",
                                           diff_date > 0 ~ "post"),
                          type = factor(type, levels = c("pre", "assess", "post"))),
                 mapping = aes(x = date, y = cases_scaling*cases, color = type),
                 alpha = 0.5) + 
      scale_colour_manual(values = c("blue", #"#00A08A", # deaths past for fitting
                                     "#CA0020"), # future predictions within 2 weeks of fit
                          labels = c("Simultaneously\npredicted cases", 
                                     "Future\npredicted cases"),
                          name = "") + 
      xlab("Date") + ylab("Daily new sypmtomatic infections\n(7 days lagged)") +
      scale_y_continuous(trans = "pseudo_log",
                         breaks = c(0, 1e1, 1e2, 1e3, 1e4, 1e5),
                         sec.axis = sec_axis(~ ./cases_scaling, name = "Daily reported cases",
                                             breaks = c(0, 1e1, 1e2, 1e3, 1e4, 1e5)/cases_scaling)) +
      scale_x_date(breaks = "2 months",
                         limits = as.Date(c("2020-02-15", NA)),
                         date_labels = "%b")}

cases_traj_fig

ggpubr::ggarrange(cases_traj_fig,  nrow = 1) %>% 
  ggsave(filename = "figures/cases_assess.pdf", width = 16, height = 7*.6)


# figure S?: counterfactuals  ----
counterfactual_sims <- readRDS("output/simulations/counterfactual_sims.Rds")

c("grey30", "#D67236", "#0b775e") %>% # some colors to use for all panels
  {list(
    {counterfactual_sims[["traj"]] %>% 
        pivot_longer(D:I_new_sympt) %>% 
        mutate(plot_col = factor(ifelse(.id == "traj", "observed", scenario), 
                                 levels = c("observed", "week_delay", "iso")),
               plot_name = ifelse(name == "D", "Cumulative Deaths", "Daily Cases")) %>%
        ggplot(aes(x = date, y = value, 
                   color = plot_col,
                   group = interaction(paramset, scenario, 
                                       mif2_iter, traj_set, .id))) + 
        geom_line(alpha = 0.1) + 
        geom_vline(xintercept = as.Date("2020-03-17"), 
                   linetype = "dashed") +
        scale_y_continuous(name = "", trans = "pseudo_log",
                           breaks = c(10^{0:4})) +
        scale_x_date(labels = date_format("%b"), date_breaks = "1 month",
                     limits = c(as.Date("2020-02-01"), NA)) +
        scale_colour_manual(values = .,
                            guide = FALSE,
                            labels = c(
                              "Reality",
                              "Shelter in Place\nStarting One Week Later", 
                              "Test and Isolate\nStarting on March 17, 2020"
                              )) +
        geom_text(data = data.frame(plot_name = c("Cumulative Deaths", "Daily Cases"),
                                    plot_label = c("a", "b"),
                                    lab_y = c(175, 2750)),
                  aes(y = lab_y, label = plot_label),
                  x = as.Date("2020-02-03"), 
                  inherit.aes = FALSE,
                  size = 8, fontface = 2) +
        facet_wrap(~plot_name, nrow = 2, scales = "free_y",
                   strip.position = "left") + 
        theme(strip.placement = "outside",
              strip.text = element_text(size = 16))},
    {counterfactual_sims[["total_D"]] %>%
      pivot_wider(names_from = scenario, values_from = D) %>% 
      mutate(delay_change = week_delay - observed, 
             iso_change = iso - observed) %>% 
      select(ends_with("change")) %>% 
      pivot_longer(ends_with("change")) %>% 
        # rbind(data.frame(name = "obs", value = NA)) %>% # add an empty row so "reality"/"obs" gets included in the legend
      mutate(name = factor(gsub("_change", "", name),
                           levels = c("obs", "delay", "iso"))) %>%
      ggplot(aes(x = value, group = name, 
                 color = name, 
                 fill = name)) +
      geom_histogram(aes(y = ..density..,),
                     position = "identity",
                     colour = "black", bins = 50,
                     lwd = 0.2, alpha = 0.5) +
      geom_density(alpha = 0.6) +
      scale_fill_manual(name = "Scenario", values = .,
                        labels = c("Reality\n",
                                   "Shelter in Place\nStarting One Week Later",
                                   "Test and Isolate\nStarting on March 17, 2020"),
                        aesthetics = c("fill", "color"),
                        drop = FALSE) +
      theme(legend.title    = element_text(size = 12),
            legend.position = c(0.67, 0.8),
            legend.text     = element_text(margin = margin(b = 10, unit = "pt"), 
                                           size = 10)) + 
      ylab("Density") +
      xlab("Difference in Total Deaths") +
      ggtitle("Estimated Mortality as of April 22, 2020") +
      annotate("text", x = -90, y = 0.032, label = "c",
               size = 8, fontface = 2) +
      geom_vline(xintercept = 0)}
    )} %>% 
  arrangeGrob(grobs = ., ncol = 2) %>%
  ggsave(filename = "figures/counterfactuals.pdf", 
         width = 10, height = 6,
         units = "in")

# figure S?: infected isoltation sensetivity ----
inf_iso_sims <- readRDS("output/simulations/inf_iso_sims.Rds")

inf_iso_sims %>% 
  group_by(inf_iso_level, background_soc_dist, scenario) %>% 
  summarise_at(vars(D), ~list(quantile(.x, probs = c(0.025, 0.05,  0.5, 0.95, 0.975)))) %>% 
  unnest_wider(D, names_sep = "_") %>% 
  rename_at(vars(starts_with("D_")), function(x) gsub("%", "", x, fixed = T)) %>% 
  ungroup %>%
  {ggplot(data = filter(., scenario != "maintain"),
          aes(x = 1 - background_soc_dist, 
             y = D_50, 
             group = 1 - inf_iso_level, 
             color = as.factor(1 - inf_iso_level))) + 
      geom_hline(yintercept = filter(., scenario == "maintain") %>% 
                   select(starts_with("D_")) %>% 
                   unlist(),
                 color = "grey", 
                 linetype = c("dotted", "dashed", "solid", "dashed", "dotted")) +
      geom_point(position = position_dodge(width = 0.05), size = 2) +
      geom_linerange(aes(ymin = D_2.5, ymax = D_97.5), 
                     position = position_dodge(width = 0.05)) + 
      geom_linerange(aes(ymin = D_5, ymax = D_95), 
                     size = 1.2,
                     position = position_dodge(width = 0.05)) + 
      scale_y_continuous(trans = "pseudo_log", 
                         breaks = c(250, 500, 1000, 2500, 5000, 10000),
                         name = "Total Deaths") + 
      scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                         name = "Background Social Distancing Effectiveness") +
      scale_color_manual(values = c("lightblue", "royalblue3", "black")
                         , labels = function(x){scales::percent(as.numeric(x))}
                         , name   = "Test-and-Isolate\nEffectiveness") +
      theme(legend.text = element_text(size = 10), 
            legend.title = element_text(size = 12),
            legend.position = c(0.85, 0.7))} %>% 
ggsave(filename = "figures/inf_iso_sensitivity.pdf", 
       width = 6, height = 5,
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
  filter(paramset <= 10, mif_step > 100) %>%
  filter(mif_step > 30 & LL_max_std > -5) %>% 
  ggplot(aes(x = mif_step, 
             # y = asinh(loglik),
             y = loglik,
             group = interaction(fit_date, paramset, mif2_iter),
             color = as.factor(paramset))) +
  xlab("Filtering Iteration") +
  ylab("log-likelihood") + 
  scale_colour_discrete(name = "Parameter Set") +
  guides(color = guide_legend(override.aes = list(alpha = 1, lwd = 2) ) ) +
  theme(
      axis.text.x = element_text(size = 17)
    , axis.text.y = element_text(size = 14)
    , axis.title.x = element_text(size = 17)
    , axis.title.y = element_text(size = 17)) +
  geom_line(alpha = 0.5)  + 
  facet_wrap(~fit_date, scales = "free")} # %>% ggsave(filename = "figures/mifs_LL_traces.pdf", width = 12, height = 6)

# violins of fit parameters ----

variable_labels <- c(
  expression(paste(Baseline~Transmission~Rate~β[1]))
, expression(paste(Initial~Number~Infected-E[0]))
, expression(paste(Shelter~In~Place~Social~Distancing-σ[SIP]))
)

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
             group = fit_date)) + 
  geom_violin(lwd = 1) + 
  xlab("Last Date of Data Used in Model Fitting") +
  ylab("Parameter Estimate") +
  theme(
      axis.text.x = element_text(size = 17)
    , axis.text.y = element_text(size = 14)
    , axis.title.x = element_text(size = 17)
    , axis.title.y = element_text(size = 17)) +
  facet_wrap(~name, scales = "free")} # %>% ggsave("figures/fit_param_violins.pdf", plot =., width = 12, height = 4)

# mif replicates traces ----
fig4_colors <- c("dodgerblue4", "firebrick4", "darkgoldenrod4")
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
  geom_line(alpha = 1, alpha = 0.8) + 
  scale_x_continuous(breaks = c(50, 150, 250)) +
  xlab("Filtering Iteration") +
  ylab("Parameter Estimate") + 
  scale_colour_manual(name = "Fitting Date", values = fig4_colors) +
  guides(color = guide_legend(override.aes = list(alpha = 1, lwd = 2) ) ) +
  theme(
      axis.text.x = element_text(size = 17)
    , axis.text.y = element_text(size = 14)
    , axis.title.x = element_text(size = 17)
    , axis.title.y = element_text(size = 17)
    , strip.text.x = element_text(size = 15)
    , strip.text.y = element_text(size = 15)) +
  facet_grid(name ~ paramset, scales = "free_y")} # %>% ggsave("figures/mif_replicate_traces.pdf", plot = ., width = 12, height = 6)

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

# Lightswitch -----

intervention_sims <- readRDS("figure_rds/lightswitch.Rds")

gg1h <- ggplot(intervention_sims %>% filter(scenario == "light", .id == "traj"), aes(date, D)) + 
  geom_line(aes(group = interaction(.id, traj_set))) +
  geom_line(
    data = intervention_sims %>% filter(scenario == "light", .id != "traj")
    , aes(group = interaction(.id, traj_set)), alpha = 0.10, colour = "grey40") +
  xlab("Date") + ylab("Cumulative Deaths")

gg2h <- ggplot(intervention_sims %>% filter(scenario == "light", .id == "traj"), aes(date, H)) + 
  geom_line(aes(group = interaction(.id, traj_set))) +
  geom_line(
    data = intervention_sims %>% filter(scenario == "light", .id != "traj")
    , aes(group = interaction(.id, traj_set)), alpha = 0.10, colour = "grey40") +
  xlab("Date") + ylab("Hospitalized")

gridExtra::grid.arrange(gg1h, gg2h, ncol = 1)

  
# Interventions -----

intervention_sims <- readRDS("figure_rds/interventions.Rds")
intervention_sims <- intervention_sims %>% mutate(total_inf = I + H)
intervention_sims$scenario <- factor(intervention_sims$scenario, levels = c("lift", "maintain", "isolate"))
fig4_colors       <- c("firebrick4", "dodgerblue4", "darkgoldenrod4")
deaths            <- deaths %>% 
  mutate(deaths_cum = deaths) %>% 
  mutate(deaths = deaths_cum - lag(deaths_cum)) %>% 
  replace_na(list(deaths = 0))

ggI1 <- ggplot(intervention_sims) + 
  geom_line(data = (intervention_sims %>% filter(.id != "traj"))
            , aes(x        = date
                  , y      = D
                  , group  = interaction(.id, scenario, traj_set)
                  , colour = scenario
            ), alpha = 0.05) + 
  geom_line(data = (intervention_sims %>% filter(.id == "traj"))
            , aes(x       = date
                  , y     = D
                  , group = interaction(traj_set, scenario)
            ), color = "grey30", alpha = 0.1) + 
  scale_x_date(labels = date_format("%b"), date_breaks = "1 month") +
  scale_colour_manual(
    name = "Intervention"
  , values = fig4_colors
  , labels = c(
    "Lift"
  , "Maintain"
  , "Test and Isolate"
    )) +
  guides(color = guide_legend(override.aes = list(alpha = 1, lwd = 2))) +
  ylab("Cumulative Deaths") +
  xlab("") +
  geom_point(data = deaths, aes(date, deaths_cum), lwd = 2, color = "black") +
  scale_y_log10(breaks = c(1, 10, 100, 1000, 5000), labels = c(1, 10, 100, 1000, 5000)) +
  theme(legend.title = element_text(size = 16)
        , axis.text.x = element_blank()
        , legend.position = c(0.8, 0.3)
        , legend.key.size = unit(.55, "cm")
        , legend.text = element_text(size = 14) 
        , plot.title = element_text(size = 12)
        , axis.ticks.x = element_blank()) +
  annotate("text", x = as.Date("2020-02-01"), y = 3000, label = "a",
           size = 8, fontface = 2)

ggI2 <- ggplot(intervention_sims) + 
  geom_line(data = (intervention_sims %>% filter(.id != "traj"))
            , aes(x        = date
                  , y      = total_inf
                  , group  = interaction(.id, scenario, traj_set)
                  , colour = scenario
            ), alpha = 0.05) + 
  geom_line(data = (intervention_sims %>% filter(.id == "traj"))
            , aes(x       = date
                  , y     = total_inf
                  , group = interaction(traj_set, scenario)
            ), color = "grey30", alpha = 0.1) + 
  scale_x_date(labels = date_format("%b"), date_breaks = "3 month") +
  scale_colour_manual(
    values = fig4_colors
  , labels = c(
    "Lift"
  , "Maintain"
  , "Test and Isolate"
    )) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12)
  , legend.title = element_text(size = 12)
  , plot.title = element_text(size = 12)) +
  ylab("Concurrent Infections") +
  xlab("Date") +
  scale_y_log10(breaks = c(1, 10, 100, 1000, 10000, 100000), labels = c(1, 10, 100, 1000, 10000, 1E5)) +
  theme(legend.title = element_text(size = 12)
        , legend.text = element_text(margin = margin(b = 10, unit = "pt"), size = 10)
        , legend.position = c(0.8, 0.3)) +
  guides(colour = FALSE) +
  annotate("text", x = as.Date("2020-02-01"), y = 200000, label = "b",
           size = 8, fontface = 2)

gridExtra::grid.arrange(ggI1, ggI2, ncol = 1)

  