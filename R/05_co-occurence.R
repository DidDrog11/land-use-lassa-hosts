# ---
# title: "05: Species Co-occurrence Analysis"
# ---

# This script calculates and visualizes the pairwise co-occurrence patterns
# among the focal species. It uses the posterior median occupancy probabilities
# from the final model to calculate Spearman's rank correlation coefficients
# for each pair of species, stratified by land use type. This produces Figure 4.

# 0. SETUP AND DATA LOADING -----------------------------------------------

source(here::here("R", "00_setup.R"))
# ggstatsploturl <- "https://cran.r-project.org/src/contrib/Archive/ggstatsplot/ggstatsplot_0.12.5.tar.gz"
# ggstatsexpressionurl <- "https://cran.r-project.org/src/contrib/Archive/statsExpressions/statsExpressions_1.6.1.tar.gz"
# install.packages(ggstatsexpressionurl, repos=NULL, type="source")
# install.packages(ggstatsploturl, repos=NULL, type="source")
library(ggstatsplot)

y_long <- read_rds(here("data", "processed_data", "y_long.rds"))
occ_covs <- read_rds(here("data", "processed_data", "occ_covs_df.rds")) %>%
  mutate(landuse = factor(landuse, levels = c("forest", "agriculture", "village")),
         village = factor(village, levels = c(village_order)),
         setting = factor(setting, levels = c("rural", "peri-urban")),
         group_landuse = factor(group_landuse, levels = c("forest - rural", "agriculture - rural", "agriculture - peri-urban",
                                                          "village - rural", "village - peri-urban")))
det_covs <- read_rds(here("data", "processed_data", "det_covs_sp.rds"))

# included species
sp_codes <- sort(unique(y_long$species))
N <- length(sp_codes)

final_model <- read_rds(here("data", "output", "model", "final_model.rds"))
final_ppc <- read_rds(here("data", "output", "model", "final_ppc.rds"))

# 1. PREPARE DATA FOR ANALYSIS --------------------------------------------

# Create a summary table of observed detections by land use. This is used later
# as a check to ensure we only run correlation tests for strata where both
# species were actually detected.
observed <- y_long %>%
  left_join(occ_covs) %>%
  janitor::tabyl(species, landuse)

observed_stratified <- y_long %>%
  left_join(occ_covs) %>%
  janitor::tabyl(species, setting, landuse)

# Extract the posterior median occupancy probability (psi) for each species at each site.
# This provides a single, representative value of occupancy for the correlation analysis.
psi_list <- as.data.frame.table(final_model$psi.samples) %>%
  mutate(Site = as.integer(Var3),
         Species = factor(Var2, labels = sp_codes)) %>%
  select(Site, Species, Freq) %>%
  group_by(Site, Species) %>%
  group_split()

psi_species <- lapply(psi_list, function(x) {
  
  samples <- 1
  
  sampled <- tibble(Site = rep(unique(x$Site), samples),
                    Species = rep(unique(x$Species), samples),
                    Psi = median(x$Freq),
                    Psi_mean = mean(x$Freq))
  
  return(sampled)
  
}) %>%
  bind_rows()


# 2. DEFINE CO-OCCURRENCE ANALYSIS FUNCTION --------------------------------

# This function takes a pair of species and calculates the Spearman's rank
# correlation between their median occupancy probabilities, stratified by land use.
cooccurrence_plot <- function(data = psi_species, species_1 = "mastomys_natalensis", species_2 = "rattus_rattus") {
  
  paired_df <- data %>%
    filter(Species %in% c(species_1, species_2)) %>%
    select(Site, Species, Psi) %>%
    pivot_wider(names_from = Species, values_from = Psi) %>%
    left_join(occ_covs, by = c("Site" = "site_code")) %>%
    mutate(group_landuse = factor(group_landuse, levels = c("forest - rural", "agriculture - rural", "agriculture - peri-urban", "village - rural", "village - peri-urban"),
                                  labels = c("Forest - Rural", "Agriculture - Rural", "Agriculture - Peri-urban", "Village - Rural", "Village - Peri-urban")))
  
  correlation <- paired_df %>%
    ggplot() +
    geom_point(aes(x = !! sym(as.character(species_1)), y = !! sym(as.character(species_2)), colour = group_landuse)) +
    scale_colour_manual(values = group_landuse_palette) +
    theme_bw() +
    labs(x = str_to_sentence(str_replace_all(species_1, "_", " ")),
         y = str_to_sentence(str_replace_all(species_2, "_", " ")),
         colour = "Stratified landuse")
  
  sp_1 <- paired_df %>%
    pull(species_1)
  sp_2 <- paired_df %>%
    pull(species_2)
  
  correlation_test_combined <- tibble(species_1 = species_1,
                                      species_2 = species_2,
                                      spearman_rho = cor.test(x = sp_1, y = sp_2, alternative = "two.sided", method = "spearman")$estimate,
                                      spearman_p = cor.test(x = sp_1, y = sp_2, alternative = "two.sided", method = "spearman")$p.value)
  
  if(observed %>% filter(species == species_1) %>% pull(forest) > 0 & 
     observed %>% filter(species == species_2) %>% pull(forest) > 0) {
    
    sp_1 <- paired_df %>%
      filter(str_detect(group_landuse, "Forest")) %>%
      pull(species_1)
    sp_2 <- paired_df %>%
      filter(str_detect(group_landuse, "Forest")) %>%
      pull(species_2)
    
    correlation_test_forest <- tibble(landuse = "Forest",
                                      species_1 = species_1,
                                      species_2 = species_2,
                                      spearman_rho = cor.test(x = sp_1, y = sp_2, alternative = "two.sided", method = "spearman")$estimate,
                                      spearman_p = cor.test(x = sp_1, y = sp_2, alternative = "two.sided", method = "spearman")$p.value)
  } else { 
    
    correlation_test_forest <- tibble(landuse = "Forest",
                                      species_1 = species_1,
                                      species_2 = species_2,
                                      spearman_rho = NA,
                                      spearman_p = NA)
    
  }
  
  if(observed %>% filter(species == species_1) %>% pull(agriculture) > 0 & 
     observed %>% filter(species == species_2) %>% pull(agriculture) > 0) {
    
    sp_1 <- paired_df %>%
      filter(str_detect(group_landuse, "Agriculture")) %>%
      pull(species_1)
    sp_2 <- paired_df %>%
      filter(str_detect(group_landuse, "Agriculture")) %>%
      pull(species_2)
    
    correlation_test_agriculture <- tibble(landuse = "Agriculture",
                                           species_1 = species_1,
                                           species_2 = species_2,
                                           spearman_rho = cor.test(x = sp_1, y = sp_2, alternative = "two.sided", method = "spearman")$estimate,
                                           spearman_p = cor.test(x = sp_1, y = sp_2, alternative = "two.sided", method = "spearman")$p.value)
    
    
  } else {
    
    correlation_test_agriculture <- tibble(landuse = "Agriculture",
                                           species_1 = species_1,
                                           species_2 = species_2,
                                           spearman_rho = NA,
                                           spearman_p = NA)
    
  }
  
  if(observed %>% filter(species == species_1) %>% pull(village) > 0 & 
     observed %>% filter(species == species_2) %>% pull(village) > 0) {
    
    sp_1 <- paired_df %>%
      filter(str_detect(group_landuse, "Village")) %>%
      pull(species_1)
    sp_2 <- paired_df %>%
      filter(str_detect(group_landuse, "Village")) %>%
      pull(species_2)
    
    correlation_test_village <- tibble(landuse = "Village",
                                       species_1 = species_1,
                                       species_2 = species_2,
                                       spearman_rho = cor.test(x = sp_1, y = sp_2, alternative = "two.sided", method = "spearman")$estimate,
                                       spearman_p = cor.test(x = sp_1, y = sp_2, alternative = "two.sided", method = "spearman")$p.value)
    
  } else {
    
    correlation_test_village <- tibble(landuse = "Village",
                                       species_1 = species_1,
                                       species_2 = species_2,
                                       spearman_rho = NA,
                                       spearman_p = NA)
    
  }
  
  correlation_tests <- bind_rows(correlation_test_forest,
                                 correlation_test_agriculture,
                                 correlation_test_village)
  
  return(list(correlation_plot = correlation,
              correlation_test_combined = correlation_test_combined,
              correlation_tests = correlation_tests))
  
}

# 3. RUN ANALYSIS FOR ALL SPECIES PAIRS -----------------------------------

# Define the list of species to include in the analysis.
species_list <- list()
sp_codes_2 <- c("mastomys_natalensis", "rattus_rattus", "mus_musculus", "crocidura_olivieri", "praomys_rostratus", "lophuromys_sikapusi", "mus_setulosus")

associations_to_test <- expand_grid(species_1 = sp_codes_2, species_2 = sp_codes_2) %>%
  mutate(species_1 = fct(species_1, levels = c("mastomys_natalensis", "praomys_rostratus", "rattus_rattus", "mus_musculus", "lophuromys_sikapusi", "mus_setulosus",
                         "crocidura_olivieri")),
         species_2 = fct(species_2, levels = c("mastomys_natalensis", "praomys_rostratus", "rattus_rattus", "mus_musculus", "lophuromys_sikapusi", "mus_setulosus",
                         "crocidura_olivieri"))) %>%
  mutate(key = paste0(species_1, species_2, sep = "")) %>%
  distinct(key, .keep_all = TRUE) %>%
  select(species_1, species_2) %>%
  group_by(species_1) %>%
  group_split()

# This loop iterates through every species pair, calls the cooccurrence_plot
# function, and stores the results in a list.
correlation_list <- vector(mode = "list", length = length(sp_codes_2))

for(n in 1:length(sp_codes_2)) {
  
  correlation_list[[n]] <- vector(mode = "list", length = nrow(associations_to_test[[n]]))
  
  species_list[[n]] <- associations_to_test[[n]]
  
  for(i in 1:nrow(species_list[[n]])) {
    
    sp_1 <- species_list[[n]][i,1] %>%
      pull()
    sp_2 <- species_list[[n]][i,2] %>%
      pull()
    
    correlation_list[[n]][[i]] <- cooccurrence_plot(species_1 = sp_1, species_2 = sp_2)
    
    
  }
}

species_order_plots <- species_order_plots[1:7]

# 4. GENERATE FIGURE 4 (Correlation Heatmap) ------------------------------

# Combine the list of results into a single data frame for plotting.
correlation_df <- lapply(correlation_list, function(x) {
  lapply(x, function(y) {
    
    y$correlation_tests
    
  })
}) %>%
  bind_rows() %>%
  filter(species_1 != species_2) %>%
  mutate(sig = case_when(spearman_p <= 0.0005 ~ TRUE,
                         is.na(spearman_p) ~ NA,
                         TRUE ~ FALSE),
         species_1 = fct_rev(species_1),
         landuse = factor(landuse, levels = c("Forest", "Agriculture", "Village")),
         strength = cut(spearman_rho, breaks = c(-1, -0.8, -0.6, -0.4, -0.2, -0.05, 0.05, 0.2, 0.4, 0.6, 0.8, 1),
                        labels = c("Very strong -ve", "Strong -ve", "Moderate -ve", "Weak -ve", "Very weak -ve", "No correlation", "Very weak +ve", "Weak +ve", "Moderate +ve", "Strong +ve", "Very strong +ve")),
         correlation_coef = round(spearman_rho, 2),
         pos_cor = case_when(correlation_coef < 0 ~ FALSE,
                             is.na(correlation_coef) ~ NA,
                             TRUE ~ TRUE),
         pair_key = pmap_chr(
           list(species_1, species_2, landuse), 
           ~ paste(sort(c(...)), collapse = "_")
         )) %>%
  distinct(pair_key, .keep_all = TRUE) %>%
  select(-pair_key) 

correlation_df_species_names <- correlation_df %>%
  mutate(
    species_1 = fct_relabel(species_1, ~ str_to_sentence(str_replace_all(., "_", " "))),
    species_2 = fct_relabel(species_2, ~ str_to_sentence(str_replace_all(., "_", " ")))
  )

correlation_df_for_plot <- correlation_df_species_names %>%
  mutate(
    correlation_label = if_else(sig == TRUE, paste0(correlation_coef, "*"), as.character(correlation_coef))
  )

# --- 2. Create the plot with two geom_label layers ---
correlation_plot <- correlation_df_for_plot %>%
  ggplot() +
  geom_tile(aes(x = species_2, y = species_1, fill = correlation_coef)) +
  
  # Layer 1: For significant AND ecologically meaningful results (bolded)
  geom_label(
    data = . %>% filter(sig == TRUE & abs(correlation_coef) > 0.4), 
    aes(x = species_2, y = species_1, label = correlation_label), 
    fontface = "bold"
  ) +
  
  # Layer 2: For all other results (not bolded)
  geom_label(
    data = . %>% filter(!(sig == TRUE & abs(correlation_coef) > 0.4)), 
    aes(x = species_2, y = species_1, label = correlation_label)
  ) +
  
  scale_fill_gradient2(
    low = "darkred", high = "darkblue", na.value = "grey", limits = c(-1, 1),
    breaks = c(1, 0.5, 0, -0.5, -1), labels = c("+1 Strong +ve", "", "None", "", "-1 Strong -ve")
  ) +
  scale_x_discrete(drop = FALSE, guide = guide_axis(n.dodge = 2)) +
  scale_y_discrete(drop = FALSE) +
  facet_wrap(~ landuse, ncol = 1) +
  labs(fill = "Strength of correlation",
       x = "Species",
       y = "Species") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.key.width = unit(1.5, units = "cm"),
        axis.text = element_text(size = 12, face = "italic"),
        strip.text = element_text(size = 12))

save_plot(correlation_plot, filename = here("output", "figures", "Figure_4.png"), base_height = 8, base_width = 8)
