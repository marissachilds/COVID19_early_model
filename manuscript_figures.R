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
  , "grid")

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

# Figure 1: R0 and Reff from different fit dates ----
# R0_covid calulation needed from ../COVID_pomp.R
# all we need from this is the R0 calculation, and we don't want to read in anything else, so lets just grab those lines
read_start = grep("covid_R0", readLines("COVID_pomp.R"))
read_end = grep("}", readLines("COVID_pomp.R"))
read_end = read_end[which(read_end > read_start)] %>% min
source(textConnection(readLines("COVID_pomp.R")[read_start:read_end]))


weekly_fits = alply(list.files("weekly_fits", full.names = T, pattern = "Santa Clara"),
                    1, 
                    function(file_path){
                      print(file_path)
                      fit <- readRDS(file_path)
                      full_Reff <- matrix(nrow = dim(fit$param_array)[1], 
                                          ncol = dim(fit$param_array)[2], 
                                          data = 0)
                      full_R0 <- matrix(nrow = dim(fit$param_array)[1], 
                                        ncol = dim(fit$param_array)[2], 
                                        data = 0)
                      # still using double forloop...                     
                      for (i in 1:dim(full_Reff)[1]) {
                        temp_dyn <- fit$dynamics_summary %>% filter(paramset == i)
                        for (j in 1:dim(full_Reff)[2]) {
                          full_Reff[i, j] <- with(fit$variable_params[i, ], covid_R0(
                            beta0est     = fit$param_array[i, j, 1]
                            , fixed_params = c(fit$fixed_params, unlist(fit$variable_params[i, ]))
                            , sd_strength  = fit$param_array[i, j, 2]
                            , prop_S = unlist(temp_dyn[temp_dyn$name == "S_now", 3]) /
                              (fit$fixed_params["N"] - temp_dyn[temp_dyn$name == "total_D", 3]))
                          )
                          full_R0[i, j] <- with(fit$variable_params[i, ], covid_R0(
                            beta0est     = fit$param_array[i, j, 1]
                            , fixed_params = c(fit$fixed_params, unlist(fit$variable_params[i, ]))
                            , sd_strength  = 1
                            , prop_S = 1)
                          )
                        }
                      }
                      fit$variable_params$county = gsub("./weekly_fit/", "", 
                                                        strsplit(file_path, split = "_")[[1]][1])
                      fit$fit_sip_eff = strsplit(file_path, split = "_")[[1]][2]
                      fit$days_minus = as.numeric(strsplit(file_path, split = "_")[[1]][4])
                      fit$county = gsub("./weekly_fit/", "", 
                                        strsplit(file_path, split = "_")[[1]][1])
                      fit$variable_params$fit_sip_eff = strsplit(file_path, split = "_")[[1]][2]
                      fit$variable_params$days_minus = as.numeric(strsplit(file_path, split = "_")[[1]][4])
                      fit$variable_params %<>% mutate(date = as.Date("2020-6-27") - days_minus)
                      fit$full_Reff = full_Reff
                      fit$full_R0 = full_R0
                      return(fit)
                    })


ldply(weekly_fits, function(fit){
  fit$variable_params %>% filter(log_lik > max(log_lik) - 2)
}) %>% 
  group_by(days_minus) %>% 
  summarise(n = n(), 
            LL_range = max(log_lik) - min(log_lik))

names(weekly_fits) = gsub(".rds", "",
                          sapply(strsplit(list.files("./weekly_fits", full.names = T, pattern = "Santa Clara"), 
                                          split = "/"), "[", 3))
Rt_time <- ldply(weekly_fits, function(fit){
  return(rbind(data.frame(name = "Reff", 
                          value = as.vector(fit$full_Reff),
                          county = fit$county,
                          fit_sip_eff = fit$fit_sip_eff,
                          days_minus = fit$days_minus,
                          stringsAsFactors = F),
               data.frame(name = "R0", 
                          value = as.vector(fit$full_R0),
                          county = fit$county,
                          fit_sip_eff = fit$fit_sip_eff,
                          days_minus = fit$days_minus,
                          stringsAsFactors = F)))
}) %>% mutate(date = as_date("2020-06-27") - as.numeric(days_minus)) 

# numbers for the paper
Rt_time %>% 
  group_by(name, date) %>% 
  summarise(
    n = n(),
    min = min(value),
    p_2.5 = quantile(value, probs = 0.025),
    median = median(value),
    # mean = mean(value),
    p_97.5 = quantile(value, probs = 0.975),
    max = max(value)
    ) %>% 
  View

breaks_fun <- function(x) {
  if (max(x) > 5) {
    seq(0, 10, by = 2)
  } else {
    seq(0, 4, 1)
  }
}
figure1 <- Rt_time %>% 
  mutate(plot_panel = "Reproduction number estimates") %>%
  ggplot() +
  geom_violin(aes(x = date,
                  y = value, fill = name, color = name,
                  group = interaction(date, name)),
              scale = "width", #"area", 
              draw_quantiles = 0.5,
              # color = "black", 
              width = 5,
              trim = TRUE, 
              position = "identity") + 
  facet_grid( plot_panel ~ ., scales  = "free_y", switch = "y", as.table = FALSE) +
  geom_point(data = epi_df %>%
              filter(date > as.Date("2020-03-01") &  date < as.Date("2020-06-27")) %>%
  mutate(plot_panel = "Deaths"),
            aes(x = date, y = deaths), color = "grey50") +
  geom_line(data = epi_df %>%
               filter(date > as.Date("2020-03-01") &  date < as.Date("2020-06-27")) %>% 
              mutate(plot_panel = "Deaths"),
             aes(x = date, y = deaths), color = "grey50", ) +
  geom_hline(data = data.frame(plot_panel = "Reproduction number estimates", y = 1), 
             aes(yintercept = y), lty = "dashed") + 
  geom_vline(data = data.frame(date = Rt_time %>% pull(date) %>% unique,
                               plot_panel = "Deaths"),
             aes(xintercept = date),
             lty = "dotted", color = "grey60") + 
  # geom_vline(xintercept = as.Date("2020-03-17")) +
  # geom_segment(data = Rt_time %>% filter(name == "Reff") %>% group_by(date) %>% summarise(ymin = min(value)),
  #              aes(x = date, xend = date, y = -2, yend = ymin),
  #              lty = "dotted") + 
  scale_color_manual(values = fig2_colors[1:2], name = "Estimate",
                     labels = c(expression(paste(R[0])),
                                expression(paste(R[eff]))),
                     aesthetics = c("colour", "fill")) + 
  scale_y_continuous(breaks = breaks_fun
    # sec.axis = sec_axis(trans = ~ (.+ 2)/0.2,
                                         # breaks = c(0, 3, 6, 9),
                                         # name = "Deaths")
                     ) +
  ylab("") +
  xlab("Last day of data used in fit") +
  theme(text=element_text(size=20),
        axis.title.y.right = element_text(hjust=0.85, vjust = 0.5),
        legend.text.align = 0,
        legend.position = c(0.08, 0.9),
        plot.margin = margin(t = 5.5, r = 5.5, b = 0, l = 5.5),
        strip.placement =  "outside") 
figure1

figure1 = ggplot_gtable(ggplot_build(figure1))
figure1$heights[7] = 2.2*figure1$heights[7]
grid.draw(figure1)

ggsave("figure1.pdf", plot = figure1,
       device = "pdf", width = 15, height = 6,
       units = "in")

# Plot 3: simulations with different fit dates ----
weekly_date<-as_date("2020-04-01")+0:12*7
weekly_path<-rev(list.files("weekly", full.names=TRUE, pattern="^weekly"))

gg_cases<-lapply(1:length(weekly_path),
                 function(i) {
                   load(weekly_path[i])
                   my_plot<-filter(SEIR.sim.f, date <= as_date("2020-07-06") & date >= as_date("2020-01-31")) %>%
                     # mutate(rep=cumsum(ifelse(diff(c(201, .id))==-200, 1, 0))) %>%
                     # mutate(case_pred = I_new_sympt) %>%
                     # mutate_at(vars("I_new_sympt"), ~log(.x+1)) %>%
                     select(.id, date, I_new_sympt, paramset) %>%
                     # add on the observed data
                     rbind(data.frame(I_new_sympt=epi_df$cases,
                                      # case_pred=epi_df$raw_cases,
                                      date=epi_df$date,
                                      .id=ifelse(epi_df$date<=weekly_date[i], "past", "future"),
                                      paramset=rep(0, nrow(epi_df))
                     )) %>%
                     mutate(fit_date = weekly_date[i]) 
                   
                   return(my_plot)
                 }
)
gg_deaths <-lapply(1:length(weekly_path),
                   function(i) {
                     load(weekly_path[i])
                     my_plot<-filter(SEIR.sim.f, date <= as_date("2020-07-06") & date >= as_date("2020-01-31")) %>%
                       # mutate(rep=cumsum(ifelse(diff(c(201, .id))==-200, 1, 0))) %>%
                       # mutate(pred_deaths=D_new) %>%
                       # mutate_at(vars("D_new"), ~log(.x+1)) %>%
                       select(.id, date, D_new, paramset) %>%
                       # add on the observed data
                       rbind(data.frame(D_new=epi_df$deaths,
                                        # pred_deaths=epi_df$raw_deaths,
                                        date=epi_df$date,
                                        .id=ifelse(epi_df$date<=weekly_date[i], "past", "future"),
                                        paramset=rep(0, nrow(epi_df))
                       )) %>%
                       mutate(fit_date = weekly_date[i]) 
                     
                     return(my_plot)
                   }
)

# combine cases and deaths to make plotting easier? 
plot3_data <- rbind(do.call(rbind, gg_cases)  %>% 
                      rename(var = I_new_sympt) %>%
                      group_by(fit_date, date) %>% 
                      {rbind(., 
                             summarise(., var = quantile(var, probs = 0.025)) %>% 
                               mutate(.id = "lower", paramset = 0),
                             summarise(., var = quantile(var, probs = 0.975)) %>% 
                               mutate(.id = "upper", paramset = 0))} %>% 
                      ungroup %>%
                      filter(.id %in% c("median", "past", "future", "lower", "upper")) %>% 
                      mutate(type = "Cases"),
                    do.call(rbind, gg_deaths)  %>% 
                      rename(var = D_new) %>% 
                      group_by(fit_date, date) %>% 
                      {rbind(., 
                             summarise(., var = quantile(var, probs = 0.025)) %>% 
                               mutate(.id = "lower", paramset = 0),
                             summarise(., var = quantile(var, probs = 0.975)) %>% 
                               mutate(.id = "upper", paramset = 0))} %>% 
                      filter(.id %in% c("median", "past", "future", "lower", "upper")) %>% 
                      mutate(type = "Deaths")) %>% 
  # separately identify future points that are beyond the 2 week comparison window used in figure 4
  mutate(.id = ifelse(date > fit_date + 14 & .id == "future", "future2", .id))

plot3 <- plot3_data %>%
  filter(date > as.Date("2020-02-15")) %>%
  filter(fit_date %in% weekly_date[seq(1, 10, by = 3)]) %>%
  mutate(type = factor(type, levels = c("Deaths", "Cases")),
         fit_date = factor(format(fit_date, format = "%B %e"),
                           levels = unique(format(fit_date, format = "%B %e")),
                           ordered = TRUE),
         var = ifelse(type == "Cases" & .id %in% c("median", "lower", "upper"), 
                                   var*0.10, var)) %>% 
  ggplot() + 
  # plotting CIs over all simulations for a fit date
  geom_ribbon(data = (. %>% filter(.id %in% c("lower", "upper")) %>% 
                        pivot_wider(names_from = .id, values_from = var)),
              aes(x = date, ymin = lower, ymax = upper),
              fill = "grey30", alpha = 0.3) +
  # plotting only the medians of simulation fits
  geom_line(data = (. %>% filter(.id == "median")), # transform the case predictions
            aes(x = date, y = var, group = interaction(.id, paramset)),
            alpha = 0.5, color = "grey30") +
  # plotting the data
  geom_point(data = (. %>% 
                       filter(.id %in% c("past", "future", "future2")) %>% 
                       mutate(color = case_when(
                         type == "Deaths" & .id == "past" ~ "A", #"fitting_data",
                         type == "Cases" & .id == "past" ~ "B", #"simulcast_cases",
                         .id == "future" ~ "C", #"future_compare",
                         .id == "future2" ~ "D" #"future_nocompare"
                       ))), 
             aes(x = date, y = var, colour = color)) +
  facet_grid(type ~ fit_date, scale = "free_y", switch = "y") + 
  geom_vline(data = expand.grid(fit_date = weekly_date[seq(1, 10, by = 3)],
                                type = c("Cases", "Deaths")) %>%
               mutate(date = fit_date,
                      fit_date = factor(format(fit_date, format = "%B %e"),
                                        levels = unique(format(fit_date, format = "%B %e")),
                                        ordered = TRUE)),
             aes(xintercept = date), linetype = "dashed") +
  scale_colour_manual(values = c("#222222", #"#00A08A", # deaths past for fitting
                                 "#1963b3", #"#134B86", #"#2166AC", #"#009999"#"#F98400" #"#46ACC8" #"#3B9AB2" #"#5BBCD6", #"#F98400", # cases past during fitting period "#F2AD00"
                                 "#CA0020", # future predictions within 2 weeks of fit
                                 alpha("#D6604D", 0.7) # "#F4A582" #  #"#EF8A62", #"#F4A582", #alpha("#FF0000", 0.4),  # future predictions > 2 weeks of fit
                                 
  ), 
  name = "",
  labels = c("Deaths used for\nmodel fitting", 
             "Simultaneously\npredicted cases",
             "Future predicted,\ncompared",
             "Future predicted,\nuncompared")) + 
  scale_x_date(breaks = as.Date(c("2020-03-01", "2020-05-01", "2020-07-01")), #date_breaks = "2 month", 
               date_labels = "%b") +
  ylab("") + 
  scale_y_continuous(sec.axis = sec_axis(~./.10, name = "predicted")) +
  xlab("") +
  theme(strip.placement = "outside",
        strip.text.x = element_text(size = 15, face = "bold"),
        strip.text.y = element_text(size = 15, face = "bold"), 
        legend.position = c(1, 0.80),
        legend.justification = c(0, 0.5),
        plot.margin = margin(t = 5.5, b = 5.5, l = 0, r = 80.5, unit = "pt")) 

plot3

new_strip_r_2 <- gtable::gtable_filter(ggplotGrob(plot3), "strip-l-2")$grobs[[1]]$grobs[[1]]
new_strip_r_2$children[[2]]$children[[1]]$label <- "Predicted Symptomatic Cases"
new_strip_r_2$children[[2]]$children[[1]]$rot <- 270

# remove the second axis from the deaths row
plot3_edit <- gtable::gtable_filter(ggplotGrob(plot3), "axis-r-1|ylab-r", 
                                   trim = FALSE, invert = TRUE)  %>%
  # add second axis label for predicted cases 
  gtable::gtable_add_grob(
    new_strip_r_2,
    t = ggplotGrob(plot3)$layout %>% filter(name == "strip-l-2") %>% pull(t), 
    l = ggplotGrob(plot3)$layout %>% filter(name == "ylab-r") %>% pull(l),
    name = "strip-r-2") %>% 
  # add rectangle around the second column of panels, which are the ones used for figure 2 simulations
  gtable::gtable_add_grob(rectGrob(gp = gpar(col = "#35274A", fill = NA, lwd = 4)),
                          t = ggplotGrob(plot3)$layout %>% filter(name == "strip-t-2") %>% pull(t), 
                          b = ggplotGrob(plot3)$layout %>% filter(name == "axis-b-2") %>% pull(b), 
                          r = 9, l = 9, z = Inf, 
                          clip = "off",
                          name = "rect")

# plot3_edit %>% reposition_legend(position='right')

grid.newpage()
grid.draw(plot3_edit)


ggsave("figure3.pdf", plot3_edit, width = 16, height = 7)

# Figure S4: all fits simulated vs observed ----

figS4 <- plot3_data %>%
  filter(date > as.Date("2020-02-15")) %>%
  mutate(type = factor(type, levels = c("Deaths", "Cases")),
         fit_date = factor(format(fit_date, format = "%B %e"),
                           levels = unique(format(fit_date, format = "%B %e")),
                           ordered = TRUE),
         var = ifelse(type == "Cases" & .id %in% c("median", "lower", "upper"), 
                      var*0.10, var)) %>% 
  ggplot() + 
  # plotting CIs over all simulations for a fit date
  geom_ribbon(data = (. %>% filter(.id %in% c("lower", "upper")) %>% 
                        pivot_wider(names_from = .id, values_from = var)),
              aes(x = date, ymin = lower, ymax = upper),
              fill = "grey30", alpha = 0.3) +
  # plotting only the medians of simulation fits
  geom_line(data = (. %>% filter(.id == "median")), # transform the case predictions
            aes(x = date, y = var, group = interaction(.id, paramset)),
            alpha = 0.5, color = "grey30") +
  # plotting the data
  geom_point(data = (. %>% 
                       filter(.id %in% c("past", "future", "future2")) %>% 
                       mutate(color = case_when(
                         type == "Deaths" & .id == "past" ~ "A", #"fitting_data",
                         type == "Cases" & .id == "past" ~ "B", #"simulcast_cases",
                         .id == "future" ~ "C", #"future_compare",
                         .id == "future2" ~ "D" #"future_nocompare"
                       ))), 
             aes(x = date, y = var, colour = color)) +
  # facet_grid(type ~ fit_date, scale = "free_y", switch = "y") + 
  geom_vline(data = expand.grid(fit_date = weekly_date,
                                type = c("Cases", "Deaths")) %>%
               mutate(date = fit_date,
                      fit_date = factor(format(fit_date, format = "%B %e"),
                                        levels = unique(format(fit_date, format = "%B %e")),
                                        ordered = TRUE)),
             aes(xintercept = date), linetype = "dashed") +
  scale_colour_manual(values = c("#222222", #"#00A08A", # deaths past for fitting
                                 "#1963b3", #"#134B86", #"#2166AC", #"#009999"#"#F98400" #"#46ACC8" #"#3B9AB2" #"#5BBCD6", #"#F98400", # cases past during fitting period "#F2AD00"
                                 "#CA0020", # future predictions within 2 weeks of fit
                                 alpha("#D6604D", 0.7) # "#F4A582" #  #"#EF8A62", #"#F4A582", #alpha("#FF0000", 0.4),  # future predictions > 2 weeks of fit
                                 
  ), 
  name = "",
  labels = c("Deaths used for\nmodel fitting", 
             "Simultaneously\npredicted cases",
             "Future predicted,\ncompared",
             "Future predicted,\nuncompared")) + 
  scale_x_date(breaks = as.Date(c("2020-03-01", "2020-05-01", "2020-07-01")), #date_breaks = "2 month", 
               date_labels = "%b") +
  ylab("") + 
  scale_y_continuous(sec.axis = sec_axis(~./.10, name = "predicted")) +
  xlab("") +
  theme(strip.placement = "outside",
        strip.text.x = element_text(size = 15, face = "bold"),
        strip.text.y = element_text(size = 15, face = "bold"), 
        legend.position = c(1, 0.80),
        legend.justification = c(0, 0.5),
        plot.margin = margin(t = 5.5, b = 5.5, l = 0, r = 80.5, unit = "pt")) 

figS4_1 <- figS4 +
  ggforce::facet_grid_paginate(type ~ fit_date, 
                               scale = "free_y", switch = "y",
                               ncol = 7, nrow = 2, page = 1)

figS4_2 <- figS4 + 
  theme(legend.position = "none") +
  ggforce::facet_grid_paginate(type ~ fit_date, 
                               scale = "free_y", switch = "y",
                               ncol = 7, nrow = 2, page = 2)
# gtable::gtable_show_layout(ggplotGrob(figS4_2))

new_strip_r_2_S4_1 <- gtable::gtable_filter(ggplotGrob(figS4_1), "strip-l-2")$grobs[[1]]$grobs[[1]]
new_strip_r_2_S4_1$children[[2]]$children[[1]]$label <- "Predicted Symptomatic Cases"
new_strip_r_2_S4_1$children[[2]]$children[[1]]$rot <- 270

new_strip_r_2_S4_2 <- gtable::gtable_filter(ggplotGrob(figS4_2), "strip-l-2")$grobs[[1]]$grobs[[1]]
new_strip_r_2_S4_2$children[[2]]$children[[1]]$label <- "Predicted Symptomatic Cases"
new_strip_r_2_S4_2$children[[2]]$children[[1]]$rot <- 270
new_strip_r_2_S4_2$children[[2]]$children[[1]]$vjust <- -2.1

# remove the second axis from the deaths row
figS4_1_edit <- gtable::gtable_filter(ggplotGrob(figS4_1), "axis-r-1|ylab-r", 
                                    trim = FALSE, invert = TRUE)  %>%
  # add second axis label for predicted cases 
  gtable::gtable_add_grob(
    new_strip_r_2_S4_1,
    t = ggplotGrob(figS4_1)$layout %>% filter(name == "strip-l-2") %>% pull(t), 
    l = ggplotGrob(figS4_1)$layout %>% filter(name == "ylab-r") %>% pull(l),
    name = "strip-r-2") %>% 
  # add rectangle around the second column of panels, which are the ones used for figure 2 simulations
  gtable::gtable_add_grob(rectGrob(gp = gpar(col = "#35274A", fill = NA, lwd = 4)),
                          t = ggplotGrob(figS4_1)$layout %>% filter(name == "strip-t-2") %>% pull(t), 
                          b = ggplotGrob(figS4_1)$layout %>% filter(name == "axis-b-2") %>% pull(b), 
                          r = 13, l = 13, z = Inf, 
                          clip = "off",
                          name = "rect")

# remove the second axis from the deaths row
figS4_2_edit <- gtable::gtable_filter(ggplotGrob(figS4_2), "axis-r-1|ylab-r", 
                                      trim = FALSE, invert = TRUE)  %>%
  # add second axis label for predicted cases 
  gtable::gtable_add_grob(
    new_strip_r_2_S4_2,
    t = ggplotGrob(figS4_2)$layout %>% filter(name == "strip-l-2") %>% pull(t), 
    l = 18, #ggplotGrob(figS4_2)$layout %>% filter(name == "ylab-r") %>% pull(l),
    name = "strip-r-2", 
    clip = "off") 


grid.newpage()
grid.draw(figS4_2_edit)

ggsave("figureS4_1.pdf", figS4_1_edit, width = 20, height = 7)
ggsave("figureS4_2.pdf", figS4_2_edit, width = 20, height = 7)

# Plot 4: MAE from different fit dates ----

mae_cases <- do.call(rbind, gg_cases) %>%
  # mutate(n=paste(.id, paramset, sep="-")) %>%
  dplyr::select(pred_cases = I_new_sympt, date, fit_date, .id, paramset) %>%
  mutate(pred_cases = pred_cases*0.1) %>%
  left_join(select(epi_df, 
                   obs_cases = cases, 
                   date=date)) %>%
  filter(date > fit_date & date <= fit_date+14) %>% 
  mutate(abs_error = abs(pred_cases - obs_cases),
         error = pred_cases - obs_cases,
         error_dir = sign(pred_cases - obs_cases)) %>%
  filter(! .id %in% c("median", "past", "future")) %>% 
  # unite(col = "group_var", c(fit_date, paramset, .id), remove = FALSE) %>% 
  group_by(fit_date, paramset, .id) %>%  # group by fit_date & simulation
  # group_by(my_group = paste(n, fit_date, sep=":")) %>%
  summarise(mean_ae = mean(abs_error),
            mean_e = mean(error),
            mean_sign = mean(error_dir),
            pct_over = sum(error_dir == 1)/length(error_dir),
            .groups = "drop") %>% 
  mutate(mean_e_scale = mean_e / mean_ae) %>% 
  # separate(my_group, c("n", "fit_date"), sep=":") %>%
  select(fit_date, mean_ae, mean_e, mean_sign, pct_over, mean_e_scale) %>%
  group_by(fit_date) %>%
  mutate(fit_mean_e = mean(mean_e_scale),
         fit_mean_sign = mean(mean_sign),
         fit_pct_over = mean(pct_over)) 
  # mutate(fit_date = format.Date(as_date(fit_date), "%b %d"))
  

mae_deaths <- do.call(rbind, gg_deaths) %>%
  # mutate(n=paste(.id, paramset, sep="-")) %>%
  dplyr::select(pred_deaths = D_new, date, fit_date, paramset, .id) %>%
  left_join(select(epi_df, 
                   obs_deaths = deaths, 
                   date=date)) %>%
  filter(! .id %in% c("median", "past", "future")) %>%
  filter(date > fit_date & date <= fit_date+14)  %>% 
  mutate(abs_error = abs(pred_deaths-obs_deaths),
         error = pred_deaths - obs_deaths, 
         error_dir = sign(pred_deaths - obs_deaths)) %>%
  # group_by(my_group = paste(n, fit_date, sep=":")) %>%
  group_by(fit_date, paramset, .id) %>%  # group by fit_date & simulation
  summarise(mean_ae = mean(abs_error),
            mean_e = mean(error),
            mean_sign = mean(error_dir),
            pct_over = sum(error_dir == 1)/length(error_dir),
            .groups = "drop") %>%
  mutate(mean_e_scale = mean_e / mean_ae) %>% 
  # separate(my_group, c("n", "fit_date"), sep=":") %>%
  select(fit_date, mean_ae, mean_e, mean_sign, pct_over, mean_e_scale) %>%
  group_by(fit_date) %>%
  mutate(fit_mean_e = mean(mean_e_scale),
         fit_mean_sign = mean(mean_sign),
         fit_pct_over = mean(pct_over)) 
  # mutate(fit_date = format.Date(as_date(fit_date), "%b %d"))
  

plot4 <- rbind(mae_cases %>% cbind(type = "Cases Mean Absolute Error"),
               mae_deaths %>% cbind(type = "Deaths Mean Absolute Error")) %>% 
  mutate_at(vars("fit_date"), ~format.Date(as_date(.x), "%b %d")) %>%
  mutate(type = factor(type, 
                       levels = c("Deaths Mean Absolute Error", 
                                  "Cases Mean Absolute Error")),
         fit_date = factor(fit_date, levels = format.Date(as_date("2020-04-01")+0:12*7, "%b %d"))) %>%
  ggplot(aes(x = fit_date, y = mean_ae, fill = fit_pct_over)) + 
  geom_violin(trim = TRUE, color = "grey30", 
              scale = "width",
              position = "identity", weight = 0.1) + 
  # geom_hline(yintercept = 0) +
  facet_wrap(~type, ncol = 1, strip.position = "left", scale = "free_y") + 
  xlab("Last day of data used in fit") +
  ylab(NULL) +
  scale_fill_gradient2(
    low = muted("blue"), mid = "white", high = muted("red"),
    midpoint = 0.5, guide = "colourbar", 
    name = "Percent over predicted",
    # breaks=c(-0.5,-0.2, 0, 0.2, 0.5),
    # labels=c("Underestimate", -0.2, 0, 0.2, "Overestimate")
    breaks = c(0, 0.25, 0.5, 0.75),
    labels =scales::percent 
  ) +
  theme(strip.background = element_blank(), strip.placement = "outside",
        strip.text = element_text(size=17)) 

plot4

ggsave("figure4.pdf", plot4, height = 8, width = 14)


# figure susceptible population over time ----