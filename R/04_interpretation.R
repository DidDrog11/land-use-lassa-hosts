# ---
# title: "04: Model Interpretation and Visualization"
# ---

# This script loads the final, best-fit multi-species occupancy model and
# generates key outputs for interpretation, including posterior predictive checks,
# coefficient plots, and the main results figure (Figure 3) showing the
# probability of occurrence across land use gradients, as well as supplementary
# figures showing the marginal effects of detection covariates.

# 0. SETUP AND DATA LOADING -----------------------------------------------

source(here::here("R", "00_setup.R"))

# Load the necessary data objects created in the previous scripts.
y_long <- read_rds(here("data", "processed_data", "y_long.rds"))
occ_covs <- read_rds(here("data", "processed_data", "occ_covs_df.rds")) %>%
  mutate(landuse = factor(landuse, levels = c("forest", "agriculture", "village")),
         village = factor(village, levels = c(village_order)),
         setting = factor(setting, levels = c("rural", "peri-urban")),
         group_landuse = factor(group_landuse, levels = c("forest - rural", "agriculture - rural", "agriculture - peri-urban",
                                                          "village - rural", "village - peri-urban")))
det_covs <- read_rds(here("data", "processed_data", "det_covs_sp.rds"))

# Define the list of species included in the final model.
sp_codes <- sort(unique(y_long$species))
N <- length(sp_codes)

# Load the final model object and its posterior predictive check object.
final_model <- read_rds(here("data", "output", "model", "final_model.rds"))
final_ppc <- read_rds(here("data", "output", "model", "final_ppc.rds"))

# 1. MODEL DIAGNOSTICS AND VALIDATION -------------------------------------

# Display a summary of the model's parameter estimates, including convergence diagnostics (R-hat, ESS).
summary(final_model)

# Visually inspect MCMC chain mixing for key parameter groups (optional).
# plot(final_model$beta.samples, density = FALSE)
# plot(final_model$alpha.samples, density = FALSE)

# Review the posterior predictive check results.
# A Bayesian p-value around 0.5 indicates adequate model fit, while values
# less than 0.1 or greater than 0.9 suggest poor fit.
summary(final_ppc)

# This function creates a scatterplot comparing the fit statistic (chi-squared)
# for the observed data versus the data replicated from the fitted model.
# Points falling along the 1:1 line indicate a good fit.
ppc_vis <- function(data = final_ppc) {
  
  ppc_df <- tibble(ppv = c(rep(1:length(data$fit.y[,1]), times = 2)),
                   "crocidura_olivieri" = c(data$fit.y[ , 1], data$fit.y.rep[, 1]),
                   "lophuromys_sikapusi" = c(data$fit.y[ ,2], data$fit.y.rep[, 2]),
                   "mastomys_natalensis" = c(data$fit.y[ ,3], data$fit.y.rep[, 3]),
                   "mus_musculus" = c(data$fit.y[ ,4], data$fit.y.rep[, 4]),
                   "mus_setulosus" = c(data$fit.y[ ,5], data$fit.y.rep[, 5]),
                   "praomys_rostratus" = c(data$fit.y[ ,6], data$fit.y.rep[, 6]),
                   "rattus_rattus" = c(data$fit.y[ ,7], data$fit.y.rep[, 7]),
                   fit = c(rep("True", times = length(data$fit.y[ , 1])), rep("Fitted", times = length(data$fit.y[ , 1])))) %>%
    arrange(ppv) %>%
    pivot_longer(cols = c(contains(sp_codes)), names_to = "Species",
                 values_to = "Fitted_value") %>%
    pivot_wider(names_from = fit, values_from = "Fitted_value") %>%
    mutate(discrepancy = case_when(True > Fitted ~ "True > Fit",
                                   TRUE ~ "Fit > True"))
  
  ppc_plot <- ggplot(ppc_df, aes(x = True, y = Fitted, colour = discrepancy)) +
    geom_point() +
    geom_abline() +
    facet_wrap(~ Species, scales = "free") +
    labs(colour = "Discrepancy",
         title = "Comparison of fit statistic for the observed data (true) and fitted model") +
    theme_bw()
  
  return(ppc_plot)
  
}

ppc_plot <- ppc_vis()

save_plot(here("output", "figures", "observed_fitted_comparison.png"), ppc_plot, base_height = 7, base_width = 8)

# Coefficient plots -------------------------------------------------------
# source(here("R", "extract_coefficient_function.R"))

# coefficients <- extract_coeff(object = final_model)

# Probability of occurrence -----------------------------------------------

# 2. GENERATE FIGURE 3 (Probability of Occurrence) ------------------------

# This function extracts the posterior samples for the probability of occurrence (psi),
# calculates the median for each site, and generates the main results plot.
species_plots <- function(data = final_model) {
  
  d0 <- as.data.frame.table(final_model$psi.samples) # Extract all the posterior draws for each site and species
  
  # Convert them into a list for each site and species Psi = the frequency
  d1 <- d0 %>%
    mutate(Site = as.integer(Var3),
           Species = factor(Var2, labels = sp_codes)) %>%
    select(Site, Species, Freq) %>%
    group_by(Site, Species) %>%
    group_split()
  
  d2 <- lapply(d1, function(x) {
    
    samples <- 1
    
    sampled <- tibble(Site = rep(unique(x$Site), samples),
                      Species = rep(unique(x$Species), samples),
                      Psi = median(x$Freq),
                      Psi_mean = mean(x$Freq))
    
    return(sampled)
    
  }) %>%
    bind_rows()
  
  landuse_df <- d2 %>%
    left_join(occ_covs, by = c("Site" = "site_code")) %>%
    mutate(Species = factor(str_to_sentence(str_replace_all(Species, "_", " ")),
                            levels = species_order_plots),
           village = factor(str_to_sentence(village), levels = str_to_sentence(village_order)),
           landuse = factor(str_to_sentence(landuse), levels = names(landuse_palette)),
           peri_urban = factor(str_to_sentence(setting), levels = c("Rural", "Peri-Urban")),
           peri_urban_landuse = factor(str_to_title(group_landuse), levels = c("Forest - Rural", "Agriculture - Rural", "Village - Rural",
                                                                               "Forest - Peri-Urban", "Agriculture - Peri-Urban", "Village - Peri-Urban")))
  
  landuse_median <- landuse_df %>%
    group_by(Species, village, landuse) %>%
    summarise(Psi = median(Psi, na.rm = TRUE)) %>%
    ungroup()
  
  landuse_median_urbanisation <- landuse_df %>%
    group_by(Species, village, landuse, peri_urban_landuse) %>%
    summarise(Psi = median(Psi, na.rm = TRUE)) %>%
    ungroup()
  
  if (!require("lemon")) install.packages("lemon")
  library(lemon)
  shift_legend3 <- function(p) {
    pnls <- cowplot::plot_to_gtable(p) %>% gtable::gtable_filter("panel") %>%
      with(setNames(grobs, layout$name)) %>% purrr::keep(~identical(.x,zeroGrob()))
    
    if( length(pnls) == 0 ) stop( "No empty facets in the plot" )
    
    lemon::reposition_legend( p, "center", panel=names(pnls) )
  }
  
  landuse_plot <- shift_legend3(landuse_df %>%
                                  ggplot() +
                                  geom_jitter(aes(y = Psi, x = landuse, colour = landuse), alpha = 0.2, height = 0) + 
                                  facet_wrap(~ Species, nrow = 2) +
                                  scale_fill_manual(values = landuse_palette, name = "Landuse") +
                                  scale_colour_manual(values = landuse_palette, name = "Landuse") +
                                  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
                                  new_scale_colour() +
                                  geom_point(data = landuse_median, aes(y = Psi, x = landuse), colour = "black", inherit.aes = FALSE) +
                                  geom_line(data = landuse_median %>%
                                              group_by(Species, village), aes(y = Psi, x = landuse, group = village, linetype = village), inherit.aes = FALSE, colour = "black") +
                                  scale_colour_manual(values = village_palette, name = "Village") +
                                  theme_bw() +
                                  guides(colour = "none") +
                                  labs(y = "Probability of occurrence (ψ)",
                                       x = "Land use",
                                       linetype = "Village"))
  
  landuse_urbanisation_plot <- shift_legend3(landuse_df %>%
                                               ggplot() +
                                               geom_jitter(aes(y = Psi, x = peri_urban_landuse, colour = landuse), alpha = 0.2) + 
                                               facet_wrap(~ Species, nrow = 2) +
                                               scale_colour_manual(values = landuse_palette, name = "Landuse") +
                                               guides(colour = guide_legend(override.aes = list(alpha = 1))) +
                                               new_scale_colour() +
                                               scale_x_discrete(drop = FALSE, labels = c("", "Rural", "", "", "Peri-Urban", "")) +
                                               geom_point(data = landuse_median_urbanisation, aes(y = Psi, x = peri_urban_landuse), colour = "black", inherit.aes = FALSE) +
                                               geom_line(data = landuse_median_urbanisation %>%
                                                           group_by(Species, village), aes(y = Psi, x = peri_urban_landuse, group = village, linetype = village), inherit.aes = FALSE, colour = "black") +
                                               scale_colour_manual(values = village_palette, name = "Village") +
                                               scale_linetype_manual(values = c("solid", "longdash", "dotted", "dotdash"), name = "Village") +
                                               theme_bw() +
                                               guides(colour = "none") +
                                               theme(axis.ticks.x = element_blank(),
                                                     strip.text = element_text(size = 12),
                                                     axis.text.x = element_text(size = 12),
                                                     legend.text = element_text(size = 12),
                                                     legend.title = element_text(size = 14)) +
                                               labs(y = "Probability of occurrence (ψ)",
                                                    x = "Study site setting"))
  
  return(list(species_data = landuse_df,
              landuse_plot = landuse_plot,
              landuse_urbanisation_plot = landuse_urbanisation_plot))
  
}

plots <- species_plots()


# Model checks

check_model <- left_join(plots$species_data, y_long %>%
                           rename(Site = site_code,
                                  Species = species) %>%
                           mutate(Species = factor(str_to_sentence(str_replace_all(Species, "_", " ")),
                                                   levels = species_order_plots)) %>%
                           group_by(Site, Species) %>%
                           summarise(observed = sum(count)),
                         by = c("Site", "Species")) %>%
  mutate(observed_bin = case_when(observed > 0 ~ 1,
                                  is.na(observed) ~ 0))

visualise_check <- ggplot(check_model) +
  geom_point(aes(x = Psi, y = observed_bin)) +
  facet_wrap(~ Species, ncol = 1)


## Species occurrence by landuse -------------------------------------------

## Species occurrence by landuse split on peri-urban/rural -----------------

save_plot(plot = plots$landuse_urbanisation_plot, filename = here("output", "figures", "Figure_3.png"), base_width = 8, base_height = 8)

# 3. GENERATE SUPP. FIGURES (Marginal Effects of Detection Covs) ----------

# This function predicts and plots the marginal effect of each detection covariate
# on the probability of detection (p), holding other covariates at their mean value.
scaling_pred <- occ_covs %>%
  mutate(scaled_distance_building = scale(distance_building)[,1],
         scaled_elevation = scale(elevation)[,1],
         bin_forest = case_when(landuse == "forest" ~ 1,
                                TRUE ~ 0),
         bin_village = case_when(landuse == "village" ~ 1,
                                 TRUE ~ 0),
         bin_lal = case_when(village == "lalehun" ~ 1,
                             TRUE ~ 0),
         bin_lam = case_when(village == "lambayama" ~ 1,
                             TRUE ~ 0),
         bin_sei = case_when(village == "seilama" ~ 1,
                             TRUE ~ 0),
         intercept = 1)

# Marginal effect of distance_building
marginal_occurrence <- function(data = plots$species_data) {
  
  list_data <- data %>%
    ungroup() %>%
    mutate(scaled_distance_building = scale(distance_building)[, 1],
           scaled_elevation = scale(elevation)[, 1]) %>%
    group_by(village) %>%
    group_split()
  
  distance_plots <- lapply(list_data, function(x) {
    x %>%
      ggplot() +
      geom_point(aes(x = scaled_distance_building, y = Psi)) +
      geom_smooth(aes(x = scaled_distance_building, y = Psi), method = "loess") +
      scale_x_continuous(n.breaks = 4) +
      facet_wrap(~ Species, nrow = 2) +
      labs(title = paste(str_to_sentence(unique(x$village))),
           x = "Scaled distance from building",
           y = "Probability of occurrence (ψ)") +
      theme_bw()
  })
  
  combined_distance_plot <- plot_grid(plotlist = distance_plots)
  
  elevation_plots <- lapply(list_data, function(x) {
    
    x %>%
      ggplot() +
      geom_point(aes(x = scaled_elevation, y = Psi)) +
      facet_wrap(~ Species, nrow = 2) +
      scale_colour_manual(values = landuse_palette) +
      labs(title = paste(str_to_sentence(unique(x$village))),
           x = "Scaled elevation",
           y = "Probability of occurrence (ψ)") +
      theme_bw()
    
  })
  
  combined_elevation_plot <- plot_grid(plotlist = elevation_plots)
  
  return(list(distance = combined_distance_plot,
              elevation = combined_elevation_plot))
  
}

occurrence_marginal_effects <- marginal_occurrence()

# Marginal effects detection --------------------------------------------------
marginal_detection <- function(data = final_model) {
  
  # Precipitation scaled values
  precipitation <- scale(c(det_covs$precipitation))
  min_precipitation <- min(precipitation, na.rm = TRUE)
  max_precipitation <- max(precipitation, na.rm = TRUE)
  pred_precipitation <- seq(from = min_precipitation, to = max_precipitation, length.out = 100)
  mean_precipitation <- mean(pred_precipitation, na.rm = TRUE)
  
  # Moon fraction scaled values
  moon <- c(det_covs$moon_fraction)
  min_moon <- min(moon, na.rm = TRUE)
  max_moon <- max(moon, na.rm = TRUE)
  pred_moon <- seq(from = min_moon, to = max_moon, length.out = 100)
  mean_moon <- mean(moon, na.rm = TRUE)
  
  # TN scaled values
  tn <- scale(c(det_covs$trap_nights))
  min_trap_nights <- min(tn, na.rm = TRUE)
  max_trap_nights <- max(tn, na.rm = TRUE)
  pred_trap_nights <- seq(from = min_trap_nights, to = max_trap_nights, length.out = 100)
  mean_trap_nights <- mean(pred_trap_nights, na.rm = TRUE)
  
  # Precipitation on detection
  precipitation_prediction <- cbind(1, pred_precipitation, mean_moon, mean_trap_nights)
  out_detection_precipitation <- predict(data, precipitation_prediction, type = "detection")
  
  d0 <- as.data.frame.table(out_detection_precipitation$p.0.samples) # Extract all the posterior draws for each value of precipitation and species
  
  # Convert them into a list for each site and species Psi = the frequency
  d1 <- d0 %>%
    mutate(precipitation = as.integer(Var3),
           Species = factor(Var2, labels = sp_codes)) %>%
    select(precipitation, Species, Freq) %>%
    group_by(precipitation, Species) %>%
    group_split()
  
  d2 <- lapply(d1, function(x) {
    
    mean = mean(x$Freq)
    lower = bayestestR::hdi(x$Freq)$CI_low
    upper = bayestestR::hdi(x$Freq)$CI_high
    
    tibble(precipitation = unique(x$precipitation),
           Species = unique(x$Species),
           mean = mean,
           lower_ci = lower,
           upper_ci = upper)
  }) %>%
    bind_rows()
  
  precipitation_df <- d2 %>%
    mutate(Species = rep(sp_codes, 100),
           scaled_precipitation = rep(pred_precipitation, each = N)) %>%
    select(-precipitation) %>%
    mutate(Species = factor(str_to_sentence(str_replace_all(Species, "_", " ")),
                            levels = species_order_plots),
           Precipitation = scaled_precipitation * attr(precipitation, "scaled:scale") + attr(precipitation, 'scaled:center'))
  
  precipitation_plot <- precipitation_df %>%
    ggplot() +
    geom_line(aes(x = Precipitation, y = mean)) +
    geom_ribbon(aes(x = Precipitation, ymin = lower_ci, ymax = upper_ci), alpha = 0.2) +
    theme_bw() + 
    scale_y_continuous(limits = c(0, 1)) + 
    facet_wrap(~ Species) + 
    labs(x = "Mean monthly rainfall (mm)", 
         y = paste0("Detection Probability (*p*)")) +
    theme(axis.title.y = ggtext::element_markdown(),
          strip.text = element_text(size = 12),
          axis.title = element_text(size = 14))
  
  # Moon on detection
  moon_prediction <- cbind(1, mean_precipitation, pred_moon, mean_trap_nights)
  out_detection_moon <- predict(data, moon_prediction, type = "detection")
  
  d0 <- as.data.frame.table(out_detection_moon$p.0.samples) # Extract all the posterior draws for each value of precipitation and species
  
  # Convert them into a list for each site and species Psi = the frequency
  d1 <- d0 %>%
    mutate(moon = as.integer(Var3),
           Species = factor(Var2, labels = sp_codes)) %>%
    select(moon, Species, Freq) %>%
    group_by(moon, Species) %>%
    group_split()
  
  d2 <- lapply(d1, function(x) {
    
    mean = mean(x$Freq)
    lower = bayestestR::hdi(x$Freq)$CI_low
    upper = bayestestR::hdi(x$Freq)$CI_high
    
    tibble(moon = unique(x$moon),
           Species = unique(x$Species),
           mean = mean,
           lower_ci = lower,
           upper_ci = upper)
  }) %>%
    bind_rows()
  
  moon_df <- d2 %>%
    mutate(Species = rep(sp_codes, 100),
           moon_fraction = rep(pred_moon, each = N)) %>%
    select(-moon) %>%
    mutate(Species = factor(str_to_sentence(str_replace_all(Species, "_", " ")),
                            levels = species_order_plots))
  
  moon_plot <- moon_df %>%
    ggplot() +
    geom_line(aes(x = moon_fraction, y = mean)) +
    geom_ribbon(aes(x = moon_fraction, ymin = lower_ci, ymax = upper_ci), alpha = 0.2) +
    theme_bw() + 
    scale_y_continuous(limits = c(0, 1)) + 
    facet_wrap(~ Species) + 
    labs(x = "Fraction of full moon", 
         y = paste0("Detection Probability (*p*)")) +
    theme(axis.title.y = ggtext::element_markdown(),
          strip.text = element_text(size = 12),
          axis.title = element_text(size = 14))
  
  # TN on detection
  tn_prediction <- cbind(1, mean_precipitation, mean_moon, pred_trap_nights)
  out_detection_tn <- predict(data, tn_prediction, type = "detection")
  
  d0 <- as.data.frame.table(out_detection_tn$p.0.samples) # Extract all the posterior draws for each value of precipitation and species
  
  # Convert them into a list for each site and species Psi = the frequency
  d1 <- d0 %>%
    mutate(tn = as.integer(Var3),
           Species = factor(Var2, labels = sp_codes)) %>%
    select(tn, Species, Freq) %>%
    group_by(tn, Species) %>%
    group_split()
  
  d2 <- lapply(d1, function(x) {
    
    mean = mean(x$Freq)
    lower = bayestestR::hdi(x$Freq)$CI_low
    upper = bayestestR::hdi(x$Freq)$CI_high
    
    tibble(tn = unique(x$tn),
           Species = unique(x$Species),
           mean = mean,
           lower_ci = lower,
           upper_ci = upper)
  }) %>%
    bind_rows()
  
  tn_df <- d2 %>%
    mutate(Species = rep(sp_codes, 100),
           scaled_tn = rep(pred_trap_nights, each = N)) %>%
    select(-tn) %>%
    mutate(Species = factor(str_to_sentence(str_replace_all(Species, "_", " ")),
                            levels = species_order_plots),
           trap_night = scaled_tn * attr(tn, "scaled:scale") + attr(tn, 'scaled:center'))
  
  tn_plot <- tn_df %>%
    ggplot() +
    geom_line(aes(x = trap_night, y = mean)) +
    geom_ribbon(aes(x = trap_night, ymin = lower_ci, ymax = upper_ci), alpha = 0.2) +
    theme_bw() + 
    scale_y_continuous(limits = c(0, 1)) + 
    facet_wrap(~ Species) + 
    labs(x = "Number of trap nights (TN)", 
         y = paste0("Detection Probability (*p*)")) +
    theme(axis.title.y = ggtext::element_markdown(),
          strip.text = element_text(size = 12),
          axis.title = element_text(size = 14))
  
  
  return(list(detection_precipitation = precipitation_plot,
              detection_moon = moon_plot,
              detection_trap_nights = tn_plot))
  
}

marginal_detection_plots <- marginal_detection()

save_plot(plot = marginal_detection_plots$detection_precipitation, filename = here("output", "figures", "Supplementary_Figure_6.png"), base_width = 8, base_height = 8)
save_plot(plot = marginal_detection_plots$detection_moon, filename = here("output", "figures", "Supplementary_Figure_7.png"), base_width = 8, base_height = 8)
save_plot(plot = marginal_detection_plots$detection_trap_nights, filename = here("output", "figures", "Supplementary_Figure_8.png"), base_width = 8, base_height = 8)
