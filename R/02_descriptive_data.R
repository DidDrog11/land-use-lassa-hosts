source(here::here("R", "00_setup.R"))

observed_data <- read_rds(here("data", "processed_data", "descriptive_data.rds"))

# This dataframe includes a row for each trap placed within a grid. The four nights of trapping observations are grouped as a single replicate. 
sites <- observed_data$sites_grids$select_site %>%
  bind_rows() %>%
  select(date_set, site_id = unique_site, site, landuse, site_easting, site_northing, village, visit, grid_number, trap_id = trap_uid, trap_easting, trap_northing, elevation)

# This dataframe consolidates rodent and trap data this may be most suitable for people interested in re-using the data
# We save this as data/output consolidated_rodent_trap.csv
observed_data$detections %>%
  rename(trap_id = trap_uid) %>%
  left_join(sites) %>%
  ungroup() %>%
  write_csv(here("data", "output", "consolidated_rodent_trap.csv"))

# This dataframe includes a row for each rodent individual trapped during the study. 
detections <- observed_data$detections %>%
  rename(trap_id = trap_uid) %>%
  left_join(sites %>%
              select(site_id, trap_id, date_set))

# This calculates the number of trap nights performed at each village study site
tn <- observed_data$non_processed_data$trap_data %>%
  tibble() %>%
  select(-geometry) %>%
  filter(village != "bambawo") %>%
  group_by(village, visit, grid_number, trap_number) %>%
  summarise(tn = n()) %>%
  mutate(trap_id = paste0(village, "_", visit, "_", grid_number, "_", trap_number)) %>%
  ungroup() %>%
  select(trap_id, tn)

tn_village_landuse <- sites %>%
  select(village, landuse, visit, trap_id) %>%
  left_join(tn, by  = "trap_id") %>%
  group_by(village, landuse, visit) %>%
  summarise(tn = sum(tn))

tn_village <- sites %>%
  select(village, landuse, visit, trap_id) %>%
  left_join(tn, by  = "trap_id") %>%
  group_by(village) %>%
  summarise(tn = sum(tn))

tn_season <- sites %>%
  select(date_set, village, landuse, visit, trap_id) %>%
  left_join(tn, by  = "trap_id")  %>%
  mutate(month = month(date_set)) %>%
  left_join(season) %>%
  group_by(village, landuse, visit, season) %>%
  summarise(tn = sum(tn))

site_id_tn <- sites %>%
  select(site_id, trap_id) %>%
  left_join(tn, by = "trap_id") %>%
  group_by(site_id) %>%
  summarise(tn = sum(tn))

# Allocate seasons to months
season <- tibble(month = 1:12, season = c(rep("Dry", 4), rep("Rainy", 6), rep("Dry", 2)))

grids <- list()

for(i in 1:length(observed_data$sites_grids$grids_for_plotting)) {
  
  grid_site <- names(observed_data$sites_grids$grids_for_plotting[i])
  
  grids[[i]] <- tibble(observed_data$sites_grids$grids_for_plotting[[i]]) %>%
    mutate(grid_id = grid_site)
  
}

grids <- grids %>%
  bind_rows() %>%
  st_as_sf() %>%
  rename(geometry = 1)

# Description trapping effort -------------------------------------------

number_trap_nights <- sum(tn$tn)

number_visits <- length(unique(sites$visit))

median_number_trap_visit <- tn_village_landuse %>%
  group_by(village, visit) %>%
  summarise(tn = sum(tn)) %>%
  group_by(village) %>%
  summarise(tn_median = median(tn),
            tn_min = min(tn),
            tn_max = max(tn))

# Raw area of traps -------------------------------------------------------
trap_locations <- sites %>%
  select(village, grid_number, visit, trap_id, trap_easting, trap_northing) %>%
  st_as_sf(coords = c("trap_easting", "trap_northing"), crs = SL_UTM) %>%
  group_by(village, grid_number, visit) %>%
  summarise()

trap_areas <- trap_locations %>%
  st_convex_hull()

trap_areas$area <- st_area(trap_areas$geometry)

tibble(trap_areas) %>% 
  summarise(mean(area))

# Description trapping locations ------------------------------------------

# plot of grids by village
grids_with_traps <- mapply(X = sites %>%
                             st_as_sf(coords = c("trap_easting", "trap_northing"), crs = SL_UTM) %>%
                             mutate(trap_number = str_split(trap_id, pattern = "_", simplify = TRUE)[, 4],
                                    village = factor(village, levels = village_order)) %>%
                             arrange(village, grid_number) %>%
                             group_by(village, grid_number) %>%
                             group_split(),
                           Y = grids %>%
                             mutate(village = factor(str_split(grid_id, "_", simplify = TRUE)[, 1], levels = village_order),
                                    grid = as.numeric(str_split(grid_id, "_", simplify = TRUE)[, 2])) %>%
                             arrange(village, grid) %>%
                             group_by(village, grid) %>%
                             group_split(),
                           function(X, Y) {
                             list(st_join(Y, X, st_contains, left = FALSE))
                           }) %>%
  do.call(rbind, .) %>%
  select(village = village.x, grid_number, site_id, landuse, grid_id, grid_polygon = geometry)

f <- function(i) paste(i, collapse = "_")

grids_with_traps$grp <- factor(sapply(st_equals(grids_with_traps), f))

grids_with_traps <- grids_with_traps %>%
  distinct() %>%
  left_join(site_id_tn, by = c("site_id")) %>%
  group_by(village) %>%
  group_split()

osm_bbox <- lapply(X = grids_with_traps, function(X) {X %>% group_by(village) %>% mutate(geometry = st_union(.)) %>% distinct(village, geometry)}) %>%
  do.call(rbind, .) %>%
  st_as_sf(crs = SL_UTM) %>%
  st_transform(crs = default_CRS) %>%
  group_split()

names(osm_bbox) <- village_order

if(!file.exists(here("data", "processed_data", "baiama_raster.tif"))) {
  
  remotes::install_github("poissonconsulting/poisspatial")
  library("poisspatial")
  
  bg_lalehun = get_googlemap(center = st_coordinates(st_centroid(osm_bbox$lalehun)), size = c(640, 640), scale = 2, format = "png8",
                             maptype = "satellite", zoom = 17) %>%
    ps_ggmap_to_raster() %>%
    rast(., crs = "EPSG:4269")
  
  crs(bg_lalehun) <- "EPSG:4269"
  
  bg_lalehun <- project(bg_lalehun, SL_UTM)
  
  RGB(bg_lalehun) <- c(1:3)
  
  plot(bg_lalehun)
  
  writeRaster(bg_lalehun, here("data", "processed_data", "lalehun_raster.tif"), overwrite = TRUE)
  
  bg_seilama = get_googlemap(center = st_coordinates(st_centroid(osm_bbox$seilama)), size = c(640, 640), scale = 2, format = "png8",
                             maptype = "satellite", zoom = 16) %>%
    ps_ggmap_to_raster() %>%
    rast(., crs = "EPSG:4269")
  
  crs(bg_seilama) <- "EPSG:4269"
  
  bg_seilama <- project(bg_seilama, SL_UTM)
  
  RGB(bg_seilama) <- c(1:3)
  
  plot(bg_seilama)
  
  writeRaster(bg_seilama, here("data", "processed_data", "seilama_raster.tif"), overwrite = TRUE)
  
  bg_baiama = get_googlemap(center = st_coordinates(st_centroid(osm_bbox$baiama)), size = c(640, 640), scale = 2, format = "png8",
                            maptype = "satellite", zoom = 15) %>%
    ps_ggmap_to_raster() %>%
    rast(., crs = "EPSG:4269")
  
  crs(bg_baiama) <- "EPSG:4269"
  
  bg_baiama <- project(bg_baiama, SL_UTM)
  
  RGB(bg_baiama) <- c(1:3)
  
  plot(bg_baiama)
  
  writeRaster(bg_baiama, here("data", "processed_data", "baiama_raster.tif"), overwrite = TRUE)
  
  bg_lambayama = get_googlemap(center = st_coordinates(st_centroid(osm_bbox$lambayama)), size = c(640, 640), scale = 2, format = "png8",
                               maptype = "satellite", zoom = 17) %>%
    ps_ggmap_to_raster() %>%
    rast(., crs = "EPSG:4269")
  
  crs(bg_lambayama) <- "EPSG:4269"
  
  bg_lambayama <- project(bg_lambayama, SL_UTM)
  
  RGB(bg_lambayama) <- c(1:3)
  
  plot(bg_lambayama)
  
  writeRaster(bg_lambayama, here("data", "processed_data", "lambayama_raster.tif"), overwrite = TRUE)
  
}

fig_1_df <- grids_with_traps

write_rds(fig_1_df, here("data", "processed_data", "fig_1_df.rds"))

# Description rodents trapped ---------------------------------------------

number_rodents <- length(unique(detections$rodent_id))

trap_success_rate <- round(number_rodents/number_trap_nights * 100, 1)

trap_success_rate_df <- detections %>%
  select(rodent_id, trap_id) %>%
  left_join(sites,
            by = "trap_id") %>%
  distinct() %>%
  group_by(village, landuse) %>%
  summarise(n_rodents = n()) %>%
  left_join(tn_village_landuse %>%
              group_by(village, landuse) %>%
              summarise(tn = sum(tn)), by = c("village", "landuse")) %>%
  mutate(trap_success = round(n_rodents/tn * 100, 1))

species_trapped <- detections %>%
  select(rodent_id, trap_id, clean_names) %>%
  left_join(sites,
            by = "trap_id") %>%
  distinct() %>%
  janitor::tabyl(clean_names) %>%
  arrange(desc(percent), desc(n), clean_names)

# Description species trapped ---------------------------------------------

species_trapped <- detections %>%
  group_by(clean_names, village) %>%
  summarise(n = n()) %>%
  drop_na(clean_names) %>%
  group_by(clean_names) %>%
  mutate(N = sum(n)) %>%
  arrange(-N)

species_order <- all_species_order
names(village_order) <- c("Baiama", "Lalehun", "Lambayama", "Seilama")
landuse_order <- c("village", "agriculture", "forest")
names(landuse_order) <- c("Village", "Agriculture", "Forest")

species_trap_success_rate <- detections %>%
  select(site_id, clean_names) %>%
  group_by(site_id, clean_names) %>%
  summarise(n = n()) %>%
  left_join(sites %>%
              tibble() %>%
              select(site_id, village, landuse),
            by = "site_id") %>%
  distinct() %>%
  group_by(clean_names, village, landuse) %>%
  summarise(n = sum(n)) %>%
  drop_na(clean_names) %>%
  left_join(., trap_success_rate_df %>%
              select(-trap_success)) %>%
  rowwise() %>%
  mutate(trap_success = round(n/tn * 100, 2),
         value = paste0(n, " (", trap_success, "%)"),
         clean_names = factor(str_to_sentence(str_replace_all(clean_names, "_", " ")), levels = species_order),
         landuse = factor(landuse, levels = landuse_order, labels = names(landuse_order)),
         village = factor(village, levels = village_order, labels = names(village_order))) %>%
  select(clean_names, village, landuse, value) %>%
  arrange(landuse) %>%
  pivot_wider(names_from = landuse, values_from = value, values_fill = "-") %>%
  arrange(clean_names, village)

species_trap_success_village <- detections %>%
  select(rodent_id, site_id, clean_names) %>%
  left_join(sites %>%
              tibble() %>%
              select(site_id, village, landuse),
            by = "site_id") %>%
  distinct()  %>%
  group_by(clean_names, village, landuse) %>%
  summarise(n = n()) %>%
  drop_na(clean_names) %>%
  left_join(tn_village) %>%
  group_by(clean_names, village, tn) %>%
  summarise(n = sum(n)) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(trap_success = round(n/tn * 100, 1),
         value = paste0(n, " (", trap_success, "%)"),
         clean_names = factor(str_to_sentence(str_replace_all(clean_names, "_", " ")), levels = species_order),
         village = factor(village, levels = village_order, labels = names(village_order)),
         Combined = value) %>%
  select(clean_names, village, Combined)

table_1b <- left_join(species_trap_success_village, species_trap_success_rate, by = c("clean_names", "village"))

# Table 1b as plot --------------------------------------------------------
# Logged detection rate

fig_2_df <- detections %>%
  left_join(sites, detections,
            by = c("village", "visit", "grid_number", "trap_id", "site_id", "date_set")) %>%
  group_by(landuse, village, clean_names) %>%
  summarise(n_detected = n()) %>%
  left_join(trap_success_rate_df %>%
              select(-trap_success), by = c("village", "landuse")) %>%
  rowwise() %>%
  mutate(detection_rate = n_detected/tn * 1000,
         clean_names = str_to_sentence(str_replace(clean_names, "_", " ")),
         village = factor(str_to_sentence(village), levels = c("Baiama", "Lalehun", "Lambayama", "Seilama")),
         landuse = factor(str_to_sentence(landuse), levels = c("Forest", "Agriculture", "Village")),
         n_detected = paste0("N = ", n_detected))

fig_2_df_reordered <- fig_2_df %>%
  mutate(clean_names = fct(clean_names, levels = c(species_order_plots)))

detection_rate_breaks <- c(0, 0.01, 0.1, 10, 20, 30)

fig_2 <- ggplot(data = fig_2_df_reordered,
                aes(x = landuse, y = fct_rev(clean_names), fill = detection_rate, label = n_detected)) +
  geom_tile() +
  geom_label(fill = "white") +
  scale_fill_viridis_c(trans = scales::log10_trans(), breaks = scales::breaks_log()) +
  facet_wrap(~ village, ncol = 2) +
  labs(y = "Species",
       fill = "Detection rate \nper 1,000 trap nights",
       x = "Landuse") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.key.width = unit(1, "cm"),
        strip.text = element_text(size = 16),
        axis.text = element_text(size = 12),
        axis.text.y = element_text(
          face = if_else(str_detect(fig_2_df_reordered$clean_names, "Mastomys natalensis"), "bold.italic", "italic")
        ))

save_plot(plot = fig_2, filename = here("output", "figures", "Figure_2_village_facet.png"), base_width = 9, base_height = 12)

fig_2_by_landuse <- detections %>%
  left_join(sites, detections,
            by = c("village", "visit", "grid_number", "trap_id", "site_id", "date_set")) %>%
  group_by(landuse, clean_names) %>%
  summarise(n_detected = n()) %>%
  left_join(trap_success_rate_df %>%
              select(-trap_success) %>%
              group_by(landuse) %>%
              summarise(tn = sum(tn)), by = c("landuse")) %>%
  rowwise() %>%
  mutate(detection_rate = n_detected/tn * 1000,
         clean_names = str_to_sentence(str_replace(clean_names, "_", " ")),
         clean_names = fct(clean_names, levels = c(species_order_plots)),
         axis_names = lvls_revalue(clean_names, 
                                  c(
                                    "Mastomys natalensis (N = 113)",
                                    "Praomys rostratus (N = 102)",
                                    "Mus musculus (N = 90)",
                                    "Rattus rattus (N = 88)",
                                    "Lophuromys sikapusi (N = 57)",
                                    "Mus setulosus (N = 43)",
                                    "Crocidura olivieri (N = 105)",
                                    "Crocidura buettikoferi (N = 23)",
                                    "Crocidura grandiceps (N = 15)",
                                    "Malacomys edwardsi (N = 11)",
                                    "Lemniscomys striatus (N = 11)",
                                    "Hylomyscus simus (N = 7)",
                                    "Hybomys planifrons (N = 7)",
                                    "Mastomys erythroleucus (N = 4)",
                                    "Crocidura theresae (N = 3)",
                                    "Gerbilliscus guineae (N = 2)",
                                    "Dasymys rufulus (N = 1)"
                                  )),
         landuse = factor(str_to_sentence(landuse), levels = c("Forest", "Agriculture", "Village")),
         n_detected = paste0("n = ", n_detected))

fig_2_landuse <- ggplot(data = fig_2_by_landuse, aes(x = landuse, y = fct_rev(axis_names), fill = detection_rate, label = n_detected)) +
  geom_tile() +
  geom_label(fill = "white") +
  scale_fill_viridis_c(trans = scales::log10_trans(), breaks = scales::breaks_log()) +
  labs(y = element_blank(),
       fill = "Detection rate \nper 1,000 trap nights",
       x = "Landuse") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.key.width = unit(1, "cm"),
        strip.text = element_text(size = 16),
        axis.text = element_text(size = 13),
        axis.text.y = element_text(
          face = if_else(levels(fct_rev(fig_2_by_landuse$axis_names)) == "Mastomys natalensis (N = 113)", "bold.italic", "italic"))
        )

save_plot(plot = fig_2_landuse, filename = here("output", "figures", "Figure_2_landuse_only.png"), base_width = 8, base_height = 7)

# Description trap rate by season ----------------------------------------------

season_detection <- detections %>%
  left_join(sites,
            by = c("village", "visit", "grid_number", "trap_id", "site_id", "date_set")) %>%
  mutate(month = month(date_set)) %>%
  left_join(season) %>%
  left_join(tn_season %>%
              group_by(village, season) %>%
              summarise(tn = sum(tn))) %>%
  group_by(clean_names, village, season, tn) %>%
  summarise(n_detected = n()) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(trap_rate = (n_detected/tn)*1000,
         clean_names = factor(str_to_sentence(str_replace(clean_names, "_", " ")), levels = all_species_order)) %>%
  ggplot() +
  geom_boxplot(aes(x = season, y = trap_rate, fill = season)) +
  facet_wrap(~ clean_names) +
  theme_bw() +
  labs(fill = element_blank(),
       y = "Detection rate (/1000 TN)",
       x = "Season",
       title = "A) Village level species detection by season")

season_detection_landuse <- detections %>%
  left_join(sites,
            by = c("village", "visit", "grid_number", "trap_id", "site_id", "date_set")) %>%
  mutate(month = month(date_set)) %>%
  left_join(season) %>%
  left_join(tn_season %>%
              group_by(village, season) %>%
              summarise(tn = sum(tn))) %>%
  group_by(clean_names, village, season, landuse, tn) %>%
  summarise(n_detected = n()) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(trap_rate = (n_detected/tn)*1000,
         clean_names = str_to_sentence(str_replace(clean_names, "_", " ")),
         landuse = factor(str_to_sentence(landuse), levels = c("Forest", "Agriculture", "Village")),
         clean_names = factor(str_to_sentence(str_replace(clean_names, "_", " ")), levels = all_species_order)) %>%
  ggplot() +
  geom_boxplot(aes(x = season, y = trap_rate, fill = landuse), position = position_dodge(preserve = "single")) +
  facet_wrap(~ clean_names) +
  scale_fill_manual(values = landuse_palette) +
  theme_bw() +
  labs(fill = element_blank(),
       y = "Detection rate (/1000 TN)",
       x = "Season",
       title = "B) Landuse level species detection by season")

save_plot(plot = season_detection, filename = here("output", "figures", "Supplementary_Figure_5a.png"), base_width = 9, base_height = 7)
save_plot(plot = season_detection_landuse, filename = here("output", "figures", "Supplementary_Figure_5b.png"), base_width = 8, base_height = 7)

# Species diversity -------------------------------------------------------

richness <- sites %>%
  tibble() %>%
  distinct(village, grid_number, site_id, landuse) %>%
  left_join(detections, by = c("village", "grid_number", "site_id")) %>%
  group_by(village, grid_number, site_id, landuse, clean_names) %>%
  summarise(n = n()) %>%
  ungroup()

richness$n[is.na(richness$clean_names)] <- 0
richness$n[!is.na(richness$clean_names)] <- 1

richness_village <- richness %>%
  select(village, clean_names, n) %>%
  filter(n != 0) %>%
  distinct() %>%
  group_by(village) %>%
  summarise(species_richness = n()) %>%
  mutate(village = factor(village, levels = village_order, labels = names(village_order)))

richness_landuse <- richness %>%
  select(landuse, clean_names, n) %>%
  filter(n != 0) %>%
  distinct() %>%
  group_by(landuse) %>%
  summarise(species_richness = n()) %>%
  mutate(landuse = factor(landuse, levels = landuse_order, labels = names(landuse_order)),
         village = "Combined")

richness_landuse_village <- richness %>%
  select(village, landuse, clean_names, n) %>%
  filter(n != 0) %>%
  distinct() %>%
  group_by(village, landuse) %>%
  summarise(species_richness = n()) %>%
  mutate(landuse = factor(landuse, levels = landuse_order, labels = names(landuse_order)),
         village = factor(village, levels = village_order, labels = names(village_order)))

richness_combined <- bind_rows(richness_village, richness_landuse, richness_landuse_village)

diversity <- sites %>%
  tibble() %>%
  distinct(village, grid_number, site_id, trap_id, landuse) %>%
  left_join(detections, by = c("village", "grid_number", "trap_id", "site_id")) %>%
  group_by(village, grid_number, site_id, landuse, clean_names) %>%
  summarise(n = sum(n())) %>%
  ungroup()

diversity$n[is.na(diversity$clean_names)] <- 0

diversity_village <- diversity %>%
  filter(!is.na(clean_names)) %>%
  group_by(village, clean_names) %>%
  summarise(n = sum(n)) %>%
  group_by(village) %>%
  summarise(N = sum(n),
            shannon_diversity = diversity(n, index = "shannon", MARGIN = 2)) %>%
  mutate(village = factor(village, levels = village_order, labels = names(village_order))) %>%
  arrange(village)

diversity_landuse <- diversity %>%
  filter(!is.na(clean_names)) %>%
  group_by(landuse, clean_names) %>%
  summarise(n = sum(n)) %>%
  group_by(landuse) %>%
  summarise(N = sum(n),
            shannon_diversity = diversity(n, index = "shannon", MARGIN = 2)) %>%
  mutate(landuse = factor(landuse, levels = landuse_order, labels = names(landuse_order)),
         village = "Combined") %>%
  arrange(landuse)

diversity_landuse_village <- diversity %>%
  filter(!is.na(clean_names)) %>%
  group_by(village, landuse, clean_names) %>%
  summarise(n = sum(n)) %>%
  group_by(village, landuse) %>%
  summarise(N = sum(n),
            shannon_diversity = diversity(n, index = "shannon", MARGIN = 2)) %>%
  arrange(-shannon_diversity) %>%
  mutate(landuse = factor(landuse, levels = landuse_order, labels = names(landuse_order)),
         village = factor(village, levels = village_order, labels = names(village_order))) %>%
  arrange(village, landuse)

diversity_combined <- bind_rows(diversity_village, diversity_landuse, diversity_landuse_village)

trap_village <- tn_village %>%
  mutate(village = factor(village, levels = village_order, labels = names(village_order)))

trap_landuse <- tn_village_landuse %>%
  group_by(landuse) %>%
  summarise(tn = sum(tn)) %>%
  mutate(landuse = factor(landuse, levels = landuse_order, labels = names(landuse_order)),
         village = "Combined")

trap_landuse_village <- tn_village_landuse %>%
  select(village, landuse, tn) %>%
  group_by(village, landuse) %>%
  summarise(tn = sum(tn)) %>%
  mutate(landuse = factor(landuse, levels = landuse_order, labels = names(landuse_order)),
         village = factor(village, levels = village_order, labels = names(village_order)))

trap_combined <- bind_rows(trap_village, trap_landuse, trap_landuse_village)


# Table 1 -----------------------------------------------------------------

table_1a <- richness_combined %>%
  left_join(diversity_combined, by = c("village", "landuse")) %>%
  left_join(trap_combined, by = c("village", "landuse")) %>%
  select(village, landuse, N, TN = tn, species_richness, shannon_diversity) %>%
  mutate(village = factor(village, levels = c("Combined", str_to_sentence(village_order)), labels = c("All villages", names(village_order))),
         landuse = case_when(is.na(landuse) ~ "Total",
                             TRUE ~ as.character(landuse)),
         landuse = factor(landuse, levels = c("Village",
                                              "Agriculture",
                                              "Forest",
                                              "Total")),
         shannon_diversity = round(shannon_diversity, 2),
         TN = paste0(TN, " (", round(N/TN * 100, 1), "%)")) %>%
  arrange(village, landuse) 

write_rds(table_1a, here("output", "tables", "table_1a.rds"))

# Supplementary Figure 3 ------------------------------------------------

# Species accumulation graphs for supplementary village level
baiama_accum <- richness %>%
  filter(village == "baiama") %>%
  group_by(site_id, clean_names) %>%
  summarise(n = sum(n)) %>%
  ungroup() %>%
  pivot_wider(names_from = clean_names, values_from = n, values_fill = 0) %>%
  select(-site_id, -`NA`) %>%
  specaccum(comm = ., method = "exact")

lalehun_accum <- richness %>%
  filter(village == "lalehun") %>%
  group_by(site_id, clean_names) %>%
  summarise(n = sum(n)) %>%
  ungroup() %>%
  pivot_wider(names_from = clean_names, values_from = n, values_fill = 0) %>%
  select(-site_id, -`NA`) %>%
  specaccum(comm = ., method = "exact")

lambayama_accum <- richness %>%
  filter(village == "lambayama") %>%
  group_by(site_id, clean_names) %>%
  summarise(n = sum(n)) %>%
  ungroup() %>%
  pivot_wider(names_from = clean_names, values_from = n, values_fill = 0) %>%
  select(-site_id, -`NA`) %>%
  specaccum(comm = ., method = "exact")

seilama_accum <- richness %>%
  filter(village == "seilama") %>%
  group_by(site_id, clean_names) %>%
  summarise(n = sum(n)) %>%
  ungroup() %>%
  pivot_wider(names_from = clean_names, values_from = n, values_fill = 0) %>%
  select(-site_id, -`NA`) %>%
  specaccum(comm = ., method = "exact")

combined_accumulation <- tibble(Village = c(rep("Baiama", length(baiama_accum$sites)),
                                            rep("Lalehun", length(lalehun_accum$sites)),
                                            rep("Lambayama", length(lambayama_accum$sites)),
                                            rep("Seilama", length(seilama_accum$sites))),
                                Sites = c(baiama_accum$sites, lalehun_accum$sites, lambayama_accum$sites, seilama_accum$sites),
                                Richness = c(baiama_accum$richness, lalehun_accum$richness, lambayama_accum$richness, seilama_accum$richness),
                                sd = c(baiama_accum$sd, lalehun_accum$sd, lambayama_accum$sd, seilama_accum$sd))

accumulation_plot <- ggplot(combined_accumulation) +
  geom_line(aes(x = Sites, y = Richness, colour = Village)) +
  geom_ribbon(aes(x = Sites, ymin = Richness - sd, ymax = Richness + sd, colour = Village, fill = Village), alpha = 0.2) +
  annotate("text", x = c(320, 660, 450, 620), y = c(10, 13, 6, 14), label = unique(combined_accumulation$Village)) +
  scale_colour_manual(values = village_palette) +
  scale_fill_manual(values = village_palette) +
  theme_bw() +
  theme(legend.position = "none")

# Species accumulation graphs for supplementary village level and landuse
baiama_accum_landuse <- richness %>%
  filter(village == "baiama") %>%
  group_by(site_id, clean_names, landuse) %>%
  summarise(n = sum(n)) %>%
  ungroup() %>%
  pivot_wider(names_from = clean_names, values_from = n, values_fill = 0) %>%
  select(-site_id, -`NA`) %>%
  group_split(landuse) %>%
  lapply(., function(x) x %>%
           select(-landuse) %>%
           specaccum(comm = ., method = "exact"))
names(baiama_accum_landuse) <- c("forest", "agriculture", "village")

lalehun_accum_landuse <- richness %>%
  filter(village == "lalehun") %>%
  group_by(site_id, clean_names, landuse) %>%
  summarise(n = sum(n)) %>%
  ungroup() %>%
  pivot_wider(names_from = clean_names, values_from = n, values_fill = 0) %>%
  select(-site_id, -`NA`) %>%
  group_split(landuse) %>%
  lapply(., function(x) x %>%
           select(-landuse) %>%
           specaccum(comm = ., method = "exact"))
names(lalehun_accum_landuse) <- c("forest", "agriculture", "village")

lambayama_accum_landuse <- richness %>%
  filter(village == "lambayama") %>%
  group_by(site_id, clean_names, landuse) %>%
  summarise(n = sum(n)) %>%
  ungroup() %>%
  pivot_wider(names_from = clean_names, values_from = n, values_fill = 0) %>%
  select(-site_id, -`NA`) %>%
  group_split(landuse) %>%
  lapply(., function(x) x %>%
           select(-landuse) %>%
           specaccum(comm = ., method = "exact"))
names(lambayama_accum_landuse) <- c("agriculture", "village")


seilama_accum_landuse <- richness %>%
  filter(village == "seilama") %>%
  group_by(site_id, clean_names, landuse) %>%
  summarise(n = sum(n)) %>%
  ungroup() %>%
  pivot_wider(names_from = clean_names, values_from = n, values_fill = 0) %>%
  select(-site_id, -`NA`) %>%
  group_split(landuse) %>%
  lapply(., function(x) x %>%
           select(-landuse) %>%
           specaccum(comm = ., method = "exact"))
names(seilama_accum_landuse) <- c("forest", "agriculture", "village")

village_landuse_accumulation <- bind_rows(
  tibble(Sites = c(baiama_accum_landuse$agriculture$sites,
                   baiama_accum_landuse$forest$sites,
                   baiama_accum_landuse$village$sites),
         Village = rep("Baiama", length(c(baiama_accum_landuse$agriculture$sites,
                                          baiama_accum_landuse$forest$sites,
                                          baiama_accum_landuse$village$sites))),
         Landuse = c(rep("Agriculture", length(baiama_accum_landuse$agriculture$sites)),
                     rep("Forest", length(baiama_accum_landuse$forest$sites)),
                     rep("Village", length(baiama_accum_landuse$village$sites))),
         Richness = c(baiama_accum_landuse$agriculture$richness,
                      baiama_accum_landuse$forest$richness,
                      baiama_accum_landuse$village$richness),
         sd = c(baiama_accum_landuse$agriculture$sd,
                baiama_accum_landuse$forest$sd,
                baiama_accum_landuse$village$sd)),
  tibble(Sites = c(lalehun_accum_landuse$agriculture$sites,
                   lalehun_accum_landuse$forest$sites,
                   lalehun_accum_landuse$village$sites),
         Village = rep("Lalehun", length(c(lalehun_accum_landuse$agriculture$sites,
                                           lalehun_accum_landuse$forest$sites,
                                           lalehun_accum_landuse$village$sites))),
         Landuse = c(rep("Agriculture", length(lalehun_accum_landuse$agriculture$sites)),
                     rep("Forest", length(lalehun_accum_landuse$forest$sites)),
                     rep("Village", length(lalehun_accum_landuse$village$sites))),
         Richness = c(lalehun_accum_landuse$agriculture$richness,
                      lalehun_accum_landuse$forest$richness,
                      lalehun_accum_landuse$village$richness),
         sd = c(lalehun_accum_landuse$agriculture$sd,
                lalehun_accum_landuse$forest$sd,
                lalehun_accum_landuse$village$sd)),
  tibble(Sites = c(lambayama_accum_landuse$agriculture$sites,
                   lambayama_accum_landuse$village$sites),
         Village = rep("Lambayama", length(c(lambayama_accum_landuse$agriculture$sites,
                                             lambayama_accum_landuse$village$sites))),
         Landuse = c(rep("Agriculture", length(lambayama_accum_landuse$agriculture$sites)),
                     rep("Village", length(lambayama_accum_landuse$village$sites))),
         Richness = c(lambayama_accum_landuse$agriculture$richness,
                      lambayama_accum_landuse$forest$richness,
                      lambayama_accum_landuse$village$richness),
         sd = c(lambayama_accum_landuse$agriculture$sd,
                lambayama_accum_landuse$forest$sd,
                lambayama_accum_landuse$village$sd)),
  tibble(Sites = c(seilama_accum_landuse$agriculture$sites,
                   seilama_accum_landuse$forest$sites,
                   seilama_accum_landuse$village$sites),
         Village = rep("Seilama", length(c(seilama_accum_landuse$agriculture$sites,
                                           seilama_accum_landuse$forest$sites,
                                           seilama_accum_landuse$village$sites))),
         Landuse = c(rep("Agriculture", length(seilama_accum_landuse$agriculture$sites)),
                     rep("Forest", length(seilama_accum_landuse$forest$sites)),
                     rep("Village", length(seilama_accum_landuse$village$sites))),
         Richness = c(seilama_accum_landuse$agriculture$richness,
                      seilama_accum_landuse$forest$richness,
                      seilama_accum_landuse$village$richness),
         sd = c(seilama_accum_landuse$agriculture$sd,
                seilama_accum_landuse$forest$sd,
                seilama_accum_landuse$village$sd)),
)

accumulation_plot_village_landuse <- village_landuse_accumulation %>%
  mutate(Landuse = factor(Landuse, levels = c("Forest", "Agriculture", "Village"))) %>%
  ggplot() +
  geom_line(aes(x = Sites, y = Richness, colour = Village)) +
  geom_ribbon(aes(x = Sites, ymin = Richness - sd, ymax = Richness + sd, colour = Village, fill = Village), alpha = 0.2) +
  facet_wrap(~ Landuse) +
  scale_colour_manual(values = village_palette) +
  scale_fill_manual(values = village_palette) +
  theme_bw()

save_plot(plot = accumulation_plot + labs(title = "A)"), base_width = 10, base_height = 8, filename = here("output", "figures", "Supplementary_Figure_3a.png"))
save_plot(plot = accumulation_plot_village_landuse  + labs(title = "B)"), base_width = 10, base_height = 8, filename = here("output", "figures", "Supplementary_Figure_3b.png"))

# Discussion description --------------------------------------------------
# Trap success in buildings, village and other for comparison
indoor_traps <- observed_data$sites_grids$select_site %>%
  bind_rows() %>%
  group_by(village, visit, trap_easting, trap_northing) %>%
  mutate(traps = n()) %>%
  filter(site_use == "housing" | (landuse == "village" & traps >= 3)) %>%
  ungroup() %>%
  mutate(trap_id = paste0(village, "_", visit, "_", grid_number, "_", trap_number)) %>%
  select(site_id = unique_site, village, visit, trap_id) %>%
  left_join(tn) %>%
  mutate(landuse = "village_inside") %>%
  group_by(village, visit, site_id, landuse) %>%
  summarise(tn = sum(tn))

indoor_traps_tn <- indoor_traps %>%
  group_by(village) %>%
  summarise(tn = sum(tn))

detections_visit_indoors <- detections %>%
  group_by(site_id, visit, village) %>%
  summarise(n = n()) %>%
  filter(site_id %in% indoor_traps$site_id) %>%
  group_by(village) %>%
  summarise(n = sum(n))

left_join(detections_visit_indoors, indoor_traps_tn) %>%
  mutate(ts = round(n/tn * 100, 1))
