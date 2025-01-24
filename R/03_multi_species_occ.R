# First step is to make the data into the required format. The load_data.R  #
# script has read in the data and performed the initial reshaping of the    #
# data. The next steps are to do this for the multi-species model           #
source(here::here("R", "00_setup.R"))

observed_data <- read_rds(here("data", "processed_data", "descriptive_data.rds"))

detections <- observed_data$detections %>%
  rename(trap_id = trap_uid,
         species = clean_names) %>%
  arrange(village, visit, grid_number)

sites <- observed_data$sites_grids$select_site %>%
  bind_rows() %>%
  mutate(village = as.character(village)) %>%
  distinct(site_id = unique_site, village, visit, grid_number, site_easting, site_northing, site, trap_id = trap_uid, date_set) %>%
  arrange(village, visit, grid_number, site)

# names of the sites trapping occurred at, 2068
site_match <- tibble(site_code = 1:length(unique(sites$site_id)),
                     site_id = unique(sites$site_id))
site_codes <- site_match$site_code

trap_nights <- read_rds(here("data", "processed_data", "trap_nights.rds"))

# Produce y --------------------------------------

# produce a long format of detections with a single record per site, visit  #
# and species                                                               #

y_long <- detections %>%
  left_join(sites %>%
              select(trap_id, site_id),
            by = "trap_id") %>%
  group_by(site_id, visit, species) %>%
  summarise(count = n())

included_species <- y_long %>%
  group_by(species) %>%
  summarise(count = sum(count)) %>%
  filter(count > 25) %>%
  pull(species)

y_long <- y_long %>%
  filter(species %in% included_species)

# names of the species trapped
sp_codes <- sort(unique(y_long$species))

# define which sites were surveyed at each replicate
trap_nights_df <- site_match %>%
  left_join(trap_nights) %>%
  arrange(visit) %>%
  mutate(trap_nights = case_when(trap_nights > 1 ~ 1,
                                 TRUE ~ as.numeric(NA))) %>%
  pivot_wider(names_from = visit, values_from = trap_nights) %>%
  arrange(site_code) %>%
  select(-site_id)

# associate sites with number in site_match
y_long <- y_long %>%
  left_join(site_match, by = "site_id") %>%
  group_by(site_code) %>%
  mutate(non_0_site_code = cur_group_id()) %>%
  arrange(site_code)

write_rds(y_long, here("data", "processed_data", "y_long.rds"))

non_0_site_code <- unique(y_long$site_code)

# number of species
N <- length(sp_codes)

# number of replicates
K <- max(sites$visit)

# number of sites
J <- length(site_codes)

if(!file.exists(here("data", "processed_data", "y_sp.rds"))) {
  
  y = array(NA, dim = c(N, J, K), dimnames = list(sp_codes, site_codes[1:J], 1:K))
  
  
  for(j in 1:J) { # loop through sites
    for(k in 1:K) { # loop through replicates
      # extract data for current site/replicate combination
      curr_df <- y_long %>%
        filter(site_code == site_codes[j],
               visit == k)
      # if plot j had a detection during replicate k curr_df will have at least 1 row, if not no rodent was observed
      if(nrow(curr_df) > 0) {
        # extract the species that were observed during this site/replicate
        curr_sp <- which(sp_codes %in% curr_df$species)
        # set value to 1 for species that were observed
        y[curr_sp, j, k] <- 1
        # and to 0 otherwise
        y[-curr_sp, j, k] <- 0
      } else 
        # if plot j was not sampled at replicate k (+1 as first column is the site code) set it to NA
        if(is.na(trap_nights_df[j, k+1])) {
          y[1:N, j, k] <- NA
        } else 
          # otherwise no rodent was trapped
        { 
          y[1:N, j, k] <- 0 
        }
    }
  }
  
  str(y)
  
  # this produces our y array which can be saved to save time on repeat runs
  
  write_rds(y, here("data", "processed_data", "y_sp.rds"))
  
} else {
  
  y <- read_rds(here("data", "processed_data", "y_sp.rds"))
  
}

# summarise the total number of observations for each species at distinct sites
observed_species <- tibble(species = names(apply(y, 1, sum, na.rm = TRUE)),
                           observed = apply(y, 1, sum, na.rm = TRUE),
                           check_number = y_long %>%
                             distinct(site_code, visit, species) %>%
                             group_by(species) %>%
                             summarise(n = n()) %>%
                             arrange(species) %>%
                             pull(n))

# Produce detection covariates --------------------------------------------
# here we add covariates that can impact the probability of detecting a rodent if it is present

if(!file.exists(here("data", "processed_data", "det_covs_sp.rds"))) {
  
  raw_det <- read_rds(here("data", "observed_data", "detection_covariates.rds")) %>%
    left_join(site_match, by = c("site_id")) %>%
    distinct(site_id, site_code, visit, precipitation, moon_fraction, trap_nights) %>%
    arrange(site_code)
  
  precip_mat <- matrix(NA, nrow = J, ncol = K)
  moon_mat <- matrix(NA, nrow = J, ncol = K)
  tn_mat <- matrix(NA, nrow = J, ncol = K)
  
  for (j in 1:J) { # Loop through sites
    for (k in 1:K) { # Loop through replicate surveys
      curr_vals <- raw_det %>%
        filter(site_code == site_codes[j], visit == k)
      # If the site was surveyed for the given replicate, 
      # extract the first date and time value. 
      if (nrow(curr_vals) > 0) {
        precip_mat[j, k] <- curr_vals$precipitation[1]
        moon_mat[j, k] <- curr_vals$moon_fraction[1] 
        tn_mat[j, k] <- curr_vals$trap_nights[1] 
      }
    } # k (replicates)
  } # j (sites) 
  
  det_covs <- list(precipitation = precip_mat,
                   moon_fraction = moon_mat,
                   trap_nights = tn_mat)
  
  write_rds(det_covs, here("data", "processed_data", "det_covs_sp.rds"))
  
} else {
  
  det_covs <- read_rds(here("data", "processed_data", "det_covs_sp.rds"))
  
}

# Produce occurrence covariates -------------------------------------------
raw_occ <- read_rds(here("data", "processed_data", "occurrence_covariates.rds")) %>%
  left_join(site_match, by = c("site_id")) %>%
  mutate(village = as.character(village),
         pop_quartile = case_when(str_detect(pop_quartile, "first") ~ 1,
                                  str_detect(pop_quartile, "second") ~ 2,
                                  str_detect(pop_quartile, "third") ~ 3,
                                  str_detect(pop_quartile, "fourth") ~ 4),
         setting = case_when(village == "lambayama" ~ "peri-urban",
                             TRUE ~ "rural"),
         landuse = as.character(landuse),
         group_landuse = case_when(village == "lambayama" & landuse == "agriculture" ~ "agriculture - peri-urban",
                                   village == "lambayama" & landuse == "village" ~ "village - peri-urban",
                                   village != "lambayama" & landuse == "agriculture" ~ "agriculture - rural",
                                   village != "lambayama" & landuse == "village" ~ "village - rural",
                                   landuse == "forest" ~ "forest - rural", 
                                   TRUE ~ landuse)) %>%
  arrange(site_code) %>%
  distinct()

write_rds(raw_occ, here("data", "processed_data", "occ_covs_df.rds"))

landuse_mat <- matrix(NA, nrow = J, ncol = 1)
village_mat <- matrix(NA, nrow = J, ncol = 1)
setting_mat <- matrix(NA, nrow = J, ncol = 1)
group_landuse_mat <- matrix(NA, nrow = J, ncol = 1)
building_mat <- matrix(NA, nrow = J, ncol = 1)
dist_village_mat <- matrix(NA, nrow = J, ncol = 1)
elevation_mat <- matrix(NA, nrow = J, ncol = 1)
population_mat <- matrix(NA, nrow = J, ncol = 1)
population_q_mat <- matrix(NA, nrow = J, ncol = 1)

for(j in 1:J) {
  landuse_mat[[j]] <- raw_occ$landuse[[j]]
  village_mat[[j]] <- raw_occ$village[[j]]
  setting_mat[[j]] <- raw_occ$setting[[j]]
  group_landuse_mat[[j]] <- raw_occ$group_landuse[[j]]
  building_mat[[j]] <- raw_occ$distance_building[[j]]
  dist_village_mat[[j]] <- raw_occ$distance_centre[[j]]
  elevation_mat[[j]] <- raw_occ$elevation[[j]]
  population_mat[[j]] <- raw_occ$population[[j]]
  population_q_mat[[j]] <- raw_occ$pop_quartile[[j]]
}

occ_covs <- list(landuse = factor(landuse_mat, levels = c("forest", "agriculture", "village")),
                 village = factor(village_mat, levels = c(village_order)),
                 setting = factor(setting_mat, levels = c("rural", "peri-urban")),
                 group_landuse = factor(group_landuse_mat, levels = c("forest - rural", "agriculture - rural", "agriculture - peri-urban",
                                                                      "village - rural", "village - peri-urban")),
                 distance_building = building_mat,
                 distance_village = dist_village_mat,
                 elevation = elevation_mat,
                 population = population_mat,
                 population_q = factor(population_q_mat, levels = c(1, 2, 3, 4), labels = c("first", "second", "third", "fourth")))

write_rds(occ_covs, here("data", "processed_data", "occ_covs_list.rds"))

# Format site coordinates -------------------------------------------------
coords <- read_rds(here("data", "processed_data", "site_coords.rds"))

coords <- tibble(X = coords[, 1],
                 Y = coords[, 2],
                 site_id = row.names(coords)) %>%
  left_join(site_match, by = c("site_id")) %>%
  arrange(site_code) %>%
  select(X, Y) %>%
  as.matrix()


# Spatial model data ------------------------------------------------------
# Distances between sites
dist_sites <- dist(coords)
# Exponential covariance model
cov_model <- "exponential"

n_factors = 1

lambda_inits <- matrix(0, N, n_factors)

diag(lambda_inits) <- 1

lambda_inits[lower.tri(lambda_inits)] <- rnorm(sum(lower.tri(lambda_inits)))

# Create list object ------------------------------------------------------

data_msom_spatial <- list(y = y,
                          occ.covs = as.data.frame(occ_covs),
                          det.covs = det_covs,
                          coords = coords)

# Multi-species occupancy model --------------------------------------------
# 
# Model outputs are large ~700mb
# Because of this they are ignored in the .gitignore file
# They also take a long time to run ~ 3h for each
# We have made the final model out_ms_spatial_2.rds available on the Open Science Framework
# It is available at this link https://osf.io/jbm6y/?view_only=8f32ac4d8659464ca468914b8f89ae99
# Place the model in data/output/model and run the code from line 420-450
# The model comparison code from line 950 should continue to work as the model comparison data has been saved

## Model formulae ----------------------------------------------------------

# Spatial model intercept only
# Occurrence
occ_ms_formula_int <- ~ 1
# Detection
det_ms_formula <- ~ scale(precipitation) + moon_fraction + scale(trap_nights)

# Spatial model 1
# Occurrence
occ_ms_formula_1 <- ~ landuse

# Spatial model 2
# Occurrence
occ_ms_formula_2 <- ~ landuse + village

# Spatial model 3
# Occurrence
occ_ms_formula_3 <- ~ group_landuse

# Spatial model 4
# Occurrence
occ_ms_formula_4 <- ~ group_landuse + scale(distance_building)

# Spatial model 5
# Occurrence
occ_ms_formula_5 <- ~ group_landuse + scale(elevation)

# Spatial model 6
# Occurrence
occ_ms_formula_6 <- ~ group_landuse + scale(distance_building) + scale(elevation)

# Spatial model 7
# Occurrence
occ_ms_formula_7 <- ~ landuse + village + scale(distance_building) + scale(elevation)

# Spatial model 8
# Occurrence
occ_ms_formula_8 <- ~ landuse + village + scale(distance_building)

# Spatial model 9
# Occurrence
occ_ms_formula_9 <- ~ landuse + village + scale(elevation)

# Initial values spatial
ms_inits_spatial <- list(alpha.comm = 0,
                         beta.comm = 0,
                         beta = 0,
                         alpha = 0,
                         tau.sq.beta = 1,
                         tau.sq.alpha = 1,
                         lambda = lambda_inits, 
                         phi = 3 / mean(dist_sites),
                         z = apply(data_msom_spatial$y, c(1, 2), max, na.rm = TRUE))

# Prior values spatial
min_dist <- min(dist_sites)
max_dist <- max(dist_sites)

ms_priors_spatial <- list(beta.comm.normal = list(mean = 0, var = 2.72),
                          alpha.comm.normal = list(mean = 0, var = 2.72), 
                          tau.sq.beta.ig = list(a = 0.1, b = 0.1), 
                          tau.sq.alpha.ig = list(a = 0.1, b = 0.1),
                          sigma.sq.ig = list(a = 2, b = 2), 
                          phi.unif = list(a = 3 / max_dist, b = 3 / min_dist))


# M0 Intercept model ---------------------------------------------------------
if(!file.exists(here("data", "output", "model", "spatial_int.rds"))) {
  
  out_ms_spatial_int <- sfMsPGOcc(occ.formula = occ_ms_formula_int, 
                                  det.formula = det_ms_formula, 
                                  data = data_msom_spatial, 
                                  inits = ms_inits_spatial, 
                                  n.batch = 1000,
                                  n.burn = 4000,
                                  batch.length = 25,
                                  n.thin = 50,
                                  accept.rate = 0.43,
                                  priors = ms_priors_spatial,
                                  n.factors = n_factors,
                                  cov.model = cov_model, 
                                  n.omp.threads = 1, 
                                  verbose = TRUE, 
                                  NNGP = TRUE,
                                  n.report = 100, 
                                  n.chains = 5,
                                  tuning = list(phi = 2))
  
  write_rds(out_ms_spatial_int, here("data", "output", "model", "spatial_int.rds"))
  gc()
  
} else {
  
  out_ms_spatial_int <- read_rds(here("data", "output", "model", "spatial_int.rds"))
  
}

summary(out_ms_spatial_int)

# Landuse model -----------------------------------------------------------

## M1 Landuse only ------------------------------------------------------------
if(!file.exists(here("data", "output", "model", "spatial_1.rds"))) {
  
  out_ms_spatial_1 <- sfMsPGOcc(occ.formula = occ_ms_formula_1, 
                                det.formula = det_ms_formula, 
                                data = data_msom_spatial, 
                                inits = ms_inits_spatial, 
                                n.batch = 1000,
                                n.burn = 4000,
                                batch.length = 25,
                                n.thin = 50,
                                accept.rate = 0.43,
                                priors = ms_priors_spatial,
                                n.factors = n_factors,
                                cov.model = cov_model, 
                                n.omp.threads = 1, 
                                verbose = TRUE, 
                                NNGP = TRUE,
                                n.report = 100, 
                                n.chains = 5,
                                tuning = list(phi = 2))
  
  write_rds(out_ms_spatial_1, here("data", "output", "model", "spatial_1.rds"))
  gc()
  
} else {
  
  out_ms_spatial_1 <- read_rds(here("data", "output", "model", "spatial_1.rds"))
  
}

summary(out_ms_spatial_1)
# Community means rhat <1.1, ESS >565
# Community variances rhat <1.1, ESS >149
# Very high values for variances

# Detection means acceptable, rhat <1.1, ESS >1000
# Detection variances acceptable, rhat <1.1, ESS > 1000

# Species means accceptable, wide confidence intervals with greatest 13, rhat generally <1.1 apart from 1.4 for praomys intercept and agriculture and crocidura agriculture
# ESS <100 for praomys

# Detection generally acceptable

# Spatial variance high rhat 2.42 with ESS of 154

## M2 Landuse and village -----------------------------------------------------
if(!file.exists(here("data", "output", "model", "spatial_2.rds"))) {
  
  out_ms_spatial_2 <- sfMsPGOcc(occ.formula = occ_ms_formula_2, 
                                det.formula = det_ms_formula, 
                                data = data_msom_spatial, 
                                inits = ms_inits_spatial, 
                                n.batch = 1000,
                                n.burn = 4000,
                                batch.length = 25,
                                n.thin = 50,
                                accept.rate = 0.43,
                                priors = ms_priors_spatial,
                                n.factors = n_factors,
                                cov.model = cov_model, 
                                n.omp.threads = 1, 
                                verbose = TRUE, 
                                NNGP = TRUE,
                                n.report = 100, 
                                n.chains = 5,
                                tuning = list(phi = 2))
  
  write_rds(out_ms_spatial_2, here("data", "output", "model", "spatial_2.rds"))
  gc()
} else {
  
  out_ms_spatial_2 <- read_rds(here("data", "output", "model", "spatial_2.rds"))
  
}

summary(out_ms_spatial_2)
# Community means rhat <1.1, ESS >659
# Community variances rhat <1.15, ESS >235
# Very high values for variances

# Detection means acceptable, rhat <1.1, ESS >494
# Detection variances acceptable, rhat <1.1, ESS > 815

# Species means acceptable, wide confidence intervals with greatest -13, rhat generally <1.1
# ESS <100 for crocidura agriculture and rattus village

# Detection generally acceptable

# Spatial variance high rhat 3.55 with ESS of 5

## M3 Landuse village type combined -------------------------------------------
if(!file.exists(here("data", "output", "model", "spatial_3.rds"))) {
  
  out_ms_spatial_3 <- sfMsPGOcc(occ.formula = occ_ms_formula_3, 
                                det.formula = det_ms_formula, 
                                data = data_msom_spatial, 
                                inits = ms_inits_spatial, 
                                n.batch = 1000,
                                n.burn = 4000,
                                batch.length = 25,
                                n.thin = 50,
                                accept.rate = 0.43,
                                priors = ms_priors_spatial,
                                n.factors = n_factors,
                                cov.model = cov_model, 
                                n.omp.threads = 1, 
                                verbose = TRUE, 
                                NNGP = TRUE,
                                n.report = 100, 
                                n.chains = 5,
                                tuning = list(phi = 2))
  
  write_rds(out_ms_spatial_3, here("data", "output", "model", "spatial_3.rds"))
  gc() 
} else {
  
  out_ms_spatial_3 <- read_rds(here("data", "output", "model", "spatial_3.rds"))
  
}

summary(out_ms_spatial_3)
# Community means rhat <1.18, ESS >170
# Community variances rhat <3.7, ESS >29
# Very high values for variances

# Detection means acceptable, rhat <1.1, ESS >1037
# Detection variances acceptable, rhat <1.1, ESS >1476

# Species means accceptable, wide confidence intervals with greatest 32, rhat generally <1.1
# Although some very high values 4.26
# ESS <100 for lots of species

# Detection generally acceptable
# But high rhat and low ESS for several

# Spatial variance high rhat 2.37 with ESS of 23

## out_ms_spatial_2 seems like it may be the best currently ##
waicOcc(out_ms_spatial_int)
waicOcc(out_ms_spatial_1)
waicOcc(out_ms_spatial_2)
waicOcc(out_ms_spatial_3)
# This is supported by the WAIC ##

# Landuse village type combined plus covariates ----------------------------

## M4 Landuse village type combined + distance --------------------------------
if(!file.exists(here("data", "output", "model", "spatial_4.rds"))) {
  
  
  out_ms_spatial_4 <- sfMsPGOcc(occ.formula = occ_ms_formula_4, 
                                det.formula = det_ms_formula, 
                                data = data_msom_spatial, 
                                inits = ms_inits_spatial, 
                                n.batch = 1000,
                                n.burn = 4000,
                                batch.length = 25,
                                n.thin = 50,
                                accept.rate = 0.43,
                                priors = ms_priors_spatial,
                                n.factors = n_factors,
                                cov.model = cov_model, 
                                n.omp.threads = 1, 
                                verbose = TRUE, 
                                NNGP = TRUE,
                                n.report = 100, 
                                n.chains = 5,
                                tuning = list(phi = 2))
  
  write_rds(out_ms_spatial_4, here("data", "output", "model", "spatial_4.rds"))
  gc()
  
} else {
  
  out_ms_spatial_4 <- read_rds(here("data", "output", "model", "spatial_4.rds"))
  
}

summary(out_ms_spatial_4)
# Community means rhat <1.03, ESS >226
# Community variances rhat <1.96, ESS >59
# Very high values for variances

# Detection means acceptable, rhat <1.1, ESS >1037
# Detection variances acceptable, rhat <1.1, ESS >1476

# Species means accceptable, wide confidence intervals with greatest 26, rhat generally <1.1
# Although some very high values 2.6
# ESS <100 for several species

# Detection generally acceptable
# But high rhat and low ESS for several

# Spatial variance high rhat 1.6 with ESS of 46

## M5 Landuse village type combined + elevation -------------------------------
if(!file.exists(here("data", "output", "model", "spatial_5.rds"))) {
  
  out_ms_spatial_5 <- sfMsPGOcc(occ.formula = occ_ms_formula_5, 
                                det.formula = det_ms_formula, 
                                data = data_msom_spatial, 
                                inits = ms_inits_spatial, 
                                n.batch = 1000,
                                n.burn = 4000,
                                batch.length = 25,
                                n.thin = 50,
                                accept.rate = 0.43,
                                priors = ms_priors_spatial,
                                n.factors = n_factors,
                                cov.model = cov_model, 
                                n.omp.threads = 1, 
                                verbose = TRUE, 
                                NNGP = TRUE,
                                n.report = 100, 
                                n.chains = 5,
                                tuning = list(phi = 2))
  
  write_rds(out_ms_spatial_5, here("data", "output", "model", "spatial_5.rds"))
  gc()
  
} else {
  
  out_ms_spatial_5 <- read_rds(here("data", "output", "model", "spatial_5.rds"))
  
}

summary(out_ms_spatial_5)
# Community means rhat <1.05, ESS >306
# Community variances rhat <1.17, ESS >228
# Very high values for variances

# Detection means acceptable, rhat <1.027, ESS >614
# Detection variances acceptable, rhat <1.1, ESS >1244

# Species means accceptable, wide confidence intervals with greatest 11, rhat generally <1.1
# Although some moderately high values 1.53
# ESS <100 for few species

# Detection generally acceptable
# rhat <1.12 and low ESS for one

# Spatial variance high rhat 6.019 with ESS of 6

## M6 Landuse village type combined + distance + elevation --------------------
if(!file.exists(here("data", "output", "model", "spatial_6.rds"))) {
  
  out_ms_spatial_6 <- sfMsPGOcc(occ.formula = occ_ms_formula_6, 
                                det.formula = det_ms_formula, 
                                data = data_msom_spatial, 
                                inits = ms_inits_spatial, 
                                n.batch = 1000,
                                n.burn = 4000,
                                batch.length = 25,
                                n.thin = 50,
                                accept.rate = 0.43,
                                priors = ms_priors_spatial,
                                n.factors = n_factors,
                                cov.model = cov_model, 
                                n.omp.threads = 1, 
                                verbose = TRUE, 
                                NNGP = TRUE,
                                n.report = 100, 
                                n.chains = 5,
                                tuning = list(phi = 2))
  
  write_rds(out_ms_spatial_6, here("data", "output", "model", "spatial_6.rds"))
  gc()
  
} else {
  
  out_ms_spatial_6 <- read_rds(here("data", "output", "model", "spatial_6.rds"))
  
}

summary(out_ms_spatial_6)
# Community means rhat <1.04, ESS >330
# Community variances rhat <1.8, ESS >105
# Very high values for variances

# Detection means acceptable, rhat <1.017, ESS >558
# Detection variances acceptable, rhat <1.09, ESS >844

# Species means accceptable, wide confidence intervals with greatest 22.5, rhat generally <1.1
# Although some moderately high values 1.53
# ESS <100 for few species

# Detection generally acceptable
# rhat <1.19 and low ESS for one

# Spatial variance high rhat 1.22 with ESS of 130

## It seems like model 6 is most acceptable ##

waicOcc(out_ms_spatial_4)
waicOcc(out_ms_spatial_5)
waicOcc(out_ms_spatial_6)

## Supported by WAIC ##


# Landuse and village + covariates ----------------------------------------

## M7 Landuse and village + distance + elevation ---------------------
if(!file.exists(here("data", "output", "model", "spatial_7.rds"))) {
  
  out_ms_spatial_7 <- sfMsPGOcc(occ.formula = occ_ms_formula_7, 
                                det.formula = det_ms_formula, 
                                data = data_msom_spatial, 
                                inits = ms_inits_spatial, 
                                n.batch = 1000,
                                n.burn = 4000,
                                batch.length = 25,
                                n.thin = 50,
                                accept.rate = 0.43,
                                priors = ms_priors_spatial,
                                n.factors = n_factors,
                                cov.model = cov_model, 
                                n.omp.threads = 1, 
                                verbose = TRUE, 
                                NNGP = TRUE,
                                n.report = 100, 
                                n.chains = 5,
                                tuning = list(phi = 2))
  
  write_rds(out_ms_spatial_7, here("data", "output", "model", "spatial_7.rds"))
  gc()
  
} else {
  
  out_ms_spatial_7 <- read_rds(here("data", "output", "model", "spatial_7.rds"))
  
}

summary(out_ms_spatial_7)

# Analogous to M4
# Community means rhat <1.04, ESS >306
# Community variances rhat <1.93, ESS >70
# Very high values for variances

# Detection means acceptable, rhat <1.016, ESS >1597
# Detection variances acceptable, rhat <1.06, ESS >494

# Species means accceptable, wide confidence intervals with greatest 21, rhat generally <1.1
# Although some moderately high values 1.5
# ESS <100 for few species

# Detection generally acceptable
# rhat <1.21 and low ESS for two

# Spatial variance high rhat 1.81 with ESS of 31
# Seems generally better fit than M4
waicOcc(out_ms_spatial_4)
waicOcc(out_ms_spatial_7)
# WAIC agrees

## M8 Landuse and village + distance ---------------------
if(!file.exists(here("data", "output", "model", "spatial_8.rds"))) {
  
  out_ms_spatial_8 <- sfMsPGOcc(occ.formula = occ_ms_formula_8, 
                                det.formula = det_ms_formula, 
                                data = data_msom_spatial, 
                                inits = ms_inits_spatial, 
                                n.batch = 1000,
                                n.burn = 4000,
                                batch.length = 25,
                                n.thin = 50,
                                accept.rate = 0.43,
                                priors = ms_priors_spatial,
                                n.factors = n_factors,
                                cov.model = cov_model, 
                                n.omp.threads = 1, 
                                verbose = TRUE, 
                                NNGP = TRUE,
                                n.report = 100, 
                                n.chains = 5,
                                tuning = list(phi = 2))
  
  write_rds(out_ms_spatial_8, here("data", "output", "model", "spatial_8.rds"))
  gc()
  
} else {
  
  out_ms_spatial_8 <- read_rds(here("data", "output", "model", "spatial_8.rds"))
  
}

summary(out_ms_spatial_8)
# Analogous to M5
# Community means rhat <1.03, ESS >393
# Community variances rhat <1.25, ESS >82
# Very high values for variances

# Detection means acceptable, rhat <1.07, ESS >362
# Detection variances acceptable, rhat <1.04, ESS >431

# Species means accceptable, wide confidence intervals with greatest 13, rhat generally <1.1
# Although some moderately high values 1.3
# ESS <100 for few species

# Detection generally acceptable
# rhat <1.27 and low ESS for 4

# Spatial variance high rhat 3.3299 with ESS of 8
# Seems better than M5
waicOcc(out_ms_spatial_5)
waicOcc(out_ms_spatial_8)
# WAIC agrees

## M9 Landuse and village + elevation ---------------------
if(!file.exists(here("data", "output", "model", "spatial_9.rds"))) {
  
  out_ms_spatial_9 <- sfMsPGOcc(occ.formula = occ_ms_formula_9, 
                                det.formula = det_ms_formula, 
                                data = data_msom_spatial, 
                                inits = ms_inits_spatial, 
                                n.batch = 1000,
                                n.burn = 4000,
                                batch.length = 25,
                                n.thin = 50,
                                accept.rate = 0.43,
                                priors = ms_priors_spatial,
                                n.factors = n_factors,
                                cov.model = cov_model, 
                                n.omp.threads = 1, 
                                verbose = TRUE, 
                                NNGP = TRUE,
                                n.report = 100, 
                                n.chains = 5,
                                tuning = list(phi = 2))
  
  write_rds(out_ms_spatial_9, here("data", "output", "model", "spatial_9.rds"))
  gc()
  
} else {
  
  out_ms_spatial_9 <- read_rds(here("data", "output", "model", "spatial_9.rds"))
  
}
summary(out_ms_spatial_9)

# Analogous to M6
# Community means rhat <1.04, ESS >406
# Community variances rhat <1.45, ESS >87
# Very high values for variances

# Detection means acceptable, rhat <1.03, ESS >793
# Detection variances acceptable, rhat <1.12, ESS >833

# Species means accceptable, wide confidence intervals with greatest 19.9, rhat generally <1.1
# Although some moderately high values 1.31
# ESS <100 for few species

# Detection generally acceptable
# rhat <1.29 and low ESS for three

# Spatial variance high rhat 4.73 with ESS of 8
# Seems slightly better than M6
waicOcc(out_ms_spatial_6)
waicOcc(out_ms_spatial_9)
# WAIC is worse


# Long runs of two lowest WAIC from each group M4-M6 and M7-M9 ------------

## M10 Landuse village type combined + distance + elevation --------------------
# Based on M2
if(!file.exists(here("data", "output", "model", "spatial_10.rds"))) {
  
  out_ms_spatial_10 <- sfMsPGOcc(occ.formula = occ_ms_formula_2, 
                                 det.formula = det_ms_formula, 
                                 data = data_msom_spatial, 
                                 inits = ms_inits_spatial, 
                                 n.batch = 1500,
                                 n.burn = 10000,
                                 batch.length = 25,
                                 n.thin = 100,
                                 accept.rate = 0.43,
                                 priors = ms_priors_spatial,
                                 n.factors = n_factors,
                                 cov.model = cov_model, 
                                 n.omp.threads = 1, 
                                 verbose = TRUE, 
                                 NNGP = TRUE,
                                 n.report = 100, 
                                 n.chains = 5,
                                 tuning = list(phi = 2))
  
  write_rds(out_ms_spatial_10, here("data", "output", "model", "spatial_10.rds"))
  gc()
  
} else {
  
  out_ms_spatial_10 <- read_rds(here("data", "output", "model", "spatial_10.rds"))
  
}

summary(out_ms_spatial_10)
# Community means rhat <1.08, ESS >233
# Community variances rhat <1.95, ESS >97
# Very high values for variances

# Detection means acceptable, rhat <1.024, ESS >448
# Detection variances acceptable, rhat <1.02, ESS >1267

# Species means accceptable, wide confidence intervals with greatest 18, rhat generally <1.1
# Although some moderately high values 2.24
# ESS <100 for few species

# Detection generally acceptable
# rhat <1.4 for 2

# Spatial variance high rhat 17.4 with ESS of 3
# Not as good as M6
waicOcc(out_ms_spatial_2)
waicOcc(out_ms_spatial_10)


## M11 Landuse and village + distance + elevation -------------------------
# Currently these aren't necessarily better, perhaps try increasing so that the number of posterior samples is the same
# Expanded version of M8
if(!file.exists(here("data", "output", "model", "spatial_11.rds"))) {
  
  out_ms_spatial_11 <- sfMsPGOcc(occ.formula = occ_ms_formula_8, 
                                 det.formula = det_ms_formula, 
                                 data = data_msom_spatial, 
                                 inits = ms_inits_spatial, 
                                 n.batch = 1500,
                                 n.burn = 10000,
                                 batch.length = 25,
                                 n.thin = 100,
                                 accept.rate = 0.43,
                                 priors = ms_priors_spatial,
                                 n.factors = n_factors,
                                 cov.model = cov_model, 
                                 n.omp.threads = 1, 
                                 verbose = TRUE, 
                                 NNGP = TRUE,
                                 n.report = 100, 
                                 n.chains = 5,
                                 tuning = list(phi = 2))
  
  write_rds(out_ms_spatial_11, here("data", "output", "model", "spatial_11.rds"))
  gc()
  
} else {
  
  out_ms_spatial_11 <- read_rds(here("data", "output", "model", "spatial_11.rds"))
  
}

summary(out_ms_spatial_11)

# Community means rhat <1.04, ESS >684
# Community variances rhat <1.72, ESS >73
# Very high values for variances

# Detection means acceptable, rhat <1.01, ESS >621
# Detection variances acceptable, rhat <1.06, ESS >529

# Species means accceptable, wide confidence intervals with greatest 21, rhat generally <1.1
# Although some moderately high values 1.6
# ESS <100 for few species

# Detection generally acceptable
# rhat <1.15 and low ESS for 1

# Spatial variance high rhat 5.35 with ESS of 5
# Seems generally better fit than M4
waicOcc(out_ms_spatial_8)
waicOcc(out_ms_spatial_11)


# Model selection ---------------------------------------------------------

if(file.exists(here("data", "output", "model", "spatial_int.rds"))) {
  
  waic <- bind_rows(waicOcc(out_ms_spatial_int),
                    waicOcc(out_ms_spatial_1),
                    waicOcc(out_ms_spatial_2),
                    waicOcc(out_ms_spatial_3),
                    waicOcc(out_ms_spatial_4),
                    waicOcc(out_ms_spatial_5),
                    waicOcc(out_ms_spatial_6),
                    waicOcc(out_ms_spatial_7),
                    waicOcc(out_ms_spatial_8),
                    waicOcc(out_ms_spatial_9),
                    waicOcc(out_ms_spatial_10),
                    waicOcc(out_ms_spatial_11))
  
  model_descriptives <- tibble(model_number = 0:11,
                               formulae = as.character(c("intercept only", occ_ms_formula_1, occ_ms_formula_2,
                                                         occ_ms_formula_3, occ_ms_formula_4, occ_ms_formula_5,
                                                         occ_ms_formula_6, occ_ms_formula_7, occ_ms_formula_8,
                                                         occ_ms_formula_9, occ_ms_formula_2, occ_ms_formula_8)))
  
  model_comparison <- bind_cols(model_descriptives, waic) %>%
    arrange(WAIC)
  
  
  # elpd expected log pointwise predictive density
  # pD effective number of parameters
  write_rds(model_comparison, here("data", "output", "model", "model_comparison.rds"))
  
} else {
  
  model_comparison <- read_rds(here("data", "output", "model", "model_comparison.rds"))
  
}

model_comparison

summary(out_ms_spatial_2)

final_model <- out_ms_spatial_2

write_rds(final_model, here("data", "output", "model", "final_model.rds"))

# PPC ---------------------------------------------------------------------

final_ppc <- ppcOcc(final_model, fit.stat = "chi-squared", group = 1)
summary(final_ppc)

# Bayesian P-value = 0.5
# All species acceptable, crocidura is at the lower limit at 0.16
write_rds(final_ppc, here("data", "output", "model", "final_ppc.rds"))
