source(here::here("R", "00_setup.R"))

# Cleaned rodent and trap data is imported #

combined_data <- list()

combined_data$trap_data <- read_csv(here("data", "input", "trap_data.csv"))
combined_data$rodent_data <- read_csv(here("data", "input", "rodent_data.csv"))

# Prepare data into sites and detections ----------------------------------

detections <- combined_data$rodent_data %>%
  # remove trap night from trap_uid as this will be accounted for as a detection covariate
  mutate(trap_uid = paste0(village, "_", visit, "_", grid_number, "_", trap_number)) %>%
  select(rodent_id = rodent_uid, village, visit, grid_number, trap_number, trap_uid, clean_names) %>%
  group_by(trap_uid, clean_names) %>%
  mutate(n = n())

# Each four night trapping activity will be considered as a single replicate  #
# The exact location of a trap varied between replicates. To incorporate      #
# replicates we will re-number the trap sites based on the closes located     #
# trap in a subsequent visit. To do this we convert coordinates to projected  #
# UTM 29N for Sierra Leone this is EPSG:32629.                                #

sites <- combined_data$trap_data %>%
  select(date_set, village, trap_uid, visit, grid_number, trap_number, habitat_group, habitat, site_use, elevation, geometry) %>%
  filter(village != "bambawo") %>% # remove Bambawo as only used for one replicate
  mutate(geometry = gsub("[c()]", "", geometry)) %>%  # Remove 'c(' and ')'
  separate(geometry, into = c("longitude", "latitude"), sep = ", ", convert = TRUE) %>%
  relocate(latitude, .before = longitude) %>%
  # 43,266 TNs and 10,900 traps
  mutate(visit = as.numeric(as.character(visit)),
         grid_number = as.numeric(as.character(grid_number)),
         grid_number = case_when(grid_number == 6 ~ 7,
                                 TRUE ~ grid_number), # combining 6 and 7 as overlap spatially
         landuse = factor(case_when(str_detect(habitat_group, "agriculture") ~ "agriculture",
                                    str_detect(habitat_group, "village") ~ "village",
                                    !village %in% c("lambayama", "baiama") & str_detect(habitat_group, "forest") ~ "forest",
                                    village == "baiama" & grid_number == 2 ~ "agriculture",
                                    village == "baiama" & grid_number == 1 ~ "forest",
                                    village == "lambayama" & grid_number == 3 ~ "agriculture"), levels = c("forest", "agriculture", "village"))) %>%
  tibble(.) %>%
  distinct(village, visit, grid_number, trap_number, landuse, longitude, latitude, .keep_all = TRUE) %>%
  st_as_sf(coords = c("longitude", "latitude")) %>%
  st_set_crs(value = default_CRS) %>%
  st_transform(crs = SL_UTM) %>%
  mutate(trap_easting = st_coordinates(geometry)[, 1],
         trap_northing = st_coordinates(geometry)[, 2])

# Assign trap locations to grid cells -------------------------------------


assign_traps_to_cells <- function(all_sites = sites) {
  
  # create list where each element is a grid from a village containing all  #
  # trap locations                                                          #
  
  select_site <- all_sites %>%
    group_by(village, grid_number) %>%
    group_split()
  
  # make a 49 m^2 grid for each site, with a 1 meter buffer from the traps  #
  
  grids <- lapply(select_site, function(x) {
    st_make_grid(st_buffer(x, 1), cellsize = 7, square = TRUE)
  })
  
  # identify the cells of the grids that contain traps                      #
  
  sites_in_grid <- mapply(function(X, Y) {
    containing_grid <- st_within(X, Y)
  },
  X = select_site,
  Y = grids)
  
  # allocate these cells to the trap locations                              #
  
  select_site <- mapply(function(X, Y) {
    list(X %>%
           mutate(site = c(unlist(Y))) %>%
           tibble() %>%
           select(-geometry))
  },
  X = select_site,
  Y = sites_in_grid)
  
  # get the centroid of each of these grid cells and append to the trap     #
  # locations. Name each grid cell uniquely based on village_grid_site      #
  
  grid_coords <- lapply(grids, function(x) {
    site_easting <- st_coordinates(st_centroid(x))[, 1]
    site_northing <- st_coordinates(st_centroid(x))[, 2]
    site <- c(1:length(x))
    
    return(tibble(site_easting = site_easting,
                  site_northing = site_northing,
                  site = site))
  })
  
  select_site <- mapply(function(X, Y) {
    
    list(X %>%
           left_join(Y, by = "site") %>%
           mutate(unique_site = paste0(village, "_", grid_number, "_", site),
                  trap_uid = paste0(village, "_", visit, "_", grid_number, "_", trap_number)))
    
  },
  X = select_site,
  Y = grid_coords)
  
  # name the element in the list based on the village and grid_number     #
  
  names(select_site) <- lapply(select_site, function(x) {
    
    name <- paste0(unique(x$village), "_", unique(x$grid_number))
    
    return(name)
  })
  
  names(grid_coords) <- lapply(select_site, function(x) {
    
    name <- paste0(unique(x$village), "_", unique(x$grid_number))
    
    return(name)
  })
  
  names(grids) <- names(grid_coords)
  
  return(list(select_site = select_site,
              grid_coords = grid_coords,
              grids_for_plotting = grids))
  
}

sites_grids <- assign_traps_to_cells(sites)

write_rds(x = list(detections = detections,
                   sites_grids = sites_grids,
                   non_processed_data = combined_data),
          here("data", "processed_data", "descriptive_data.rds"))

sites_in_grid <- sites_grids$select_site
grid_coords <- sites_grids$grid_coords

for(i in 1:length(grid_coords)) {
  
  village_grid <- names(grid_coords[i])
  
  grid_coords[[i]] <- grid_coords[[i]] %>%
    mutate(unique_site = paste0(village_grid, "_", site))
  
}

grid_coords <- bind_rows(grid_coords)


# Visualising the locations of trapped cells in grids ---------------------

# we can visualise the location of the newly produced sites below       #

visualise_sites_in_grid <- lapply(sites_in_grid, function(x) {
  
  x %>% 
    group_by(village, grid_number, visit, site, site_easting, site_northing) %>%
    summarise(TN = n() * 4) %>%
    st_as_sf(coords = c("site_easting", "site_northing")) %>%
    st_set_crs(value = SL_UTM) %>%
    ggplot() +
    geom_sf(aes(colour = TN)) +
    facet_wrap(~ visit) +
    theme_bw() +
    labs(title = paste0(str_to_sentence(x$village), ": ", x$grid_number))
  
})

# Detection covariates ----------------------------------------------------

# Trap nights -------------------------------------------------------------
# Multiple traps are allocated to a single grid cell and 4 trap nights were conducted per trap
# We will use this as a measure of effort for detection

trap_nights <- lapply(sites_in_grid, function(X) {
  
  X <- X %>%
    group_by(unique_site, visit) %>%
    summarise(n_traps = n()) %>%
    mutate(trap_nights = n_traps * 4)
  
  return(X)
}) %>%
  bind_rows() %>%
  ungroup() %>%
  select(unique_site, visit, trap_nights)

write_rds(trap_nights %>%
            rename(site_id = unique_site), here("data", "processed_data", "trap_nights.rds"))

# Add date_set to the sites dataset to calculate moon and rainfall    #
# Date set for all sites                                              #
date_set <- bind_rows(sites_in_grid) %>%
  select(date_set, village, visit, site_id = unique_site, site_easting, site_northing, trap_id = trap_uid)

write_rds(date_set, here("data", "processed_data", "date_set.rds"))

# Monthly rainfall --------------------------------------------------------
# For the precipitation data we need to provide lon and lat points.     #
# We extract the centre of the trapping sites for this.                 #
spatial_traps <- bind_rows(sites_in_grid) %>%
  st_as_sf(coords = c("site_easting", "site_northing")) %>%
  st_set_crs(value = SL_UTM)

tile_coords <- as.data.frame(st_coordinates(st_transform(spatial_traps, crs = default_CRS))) %>%
  summarise(X = median(X),
            Y = median(Y))

# Download worldclim data for the trap coordinates
worldclim_tile(var = "prec", res = 0.5, lon = tile_coords[1], lat = tile_coords[2], path = here("data", "geodata"))

precip_rast <- rast(here("data", "geodata", "wc2.1_tiles", "tile_30_wc2.1_30s_prec.tif"))

month_split <- date_set %>%
  group_by(site_id, visit, site_easting, site_northing) %>%
  filter(date_set == min(date_set)) %>%
  mutate(month = month(date_set)) %>%
  st_as_sf(coords = c("site_easting", "site_northing"), crs = SL_UTM) %>%
  st_transform(crs = default_CRS) %>%
  group_by(month) %>%
  group_split()

month_split <- lapply(function(X) {
  
  # Extract month of interest
  month <- unique(X$month)
  # Subset monthly raster to month of interest
  precip_month <- precip_rast[[month]]
  # Convert spatial DF to vect for speed
  vect_X <- vect(X)
  # Append precipitation to month_split DF
  X$precipitation <- terra::extract(precip_month, vect_X)[, 2]
  
  return(X)
}, X = month_split)

# Moon phase --------------------------------------------------------------

month_split <- lapply(function(X) {
  
  date_trap <- ymd(unique(X$date_set))
  
  moon_fraction <- getMoonIllumination(date = date_trap) %>%
    select(date, fraction)
  
  X <- X %>%
    left_join(moon_fraction, by = c("date_set" = "date")) %>%
    rename(moon_fraction = fraction)
  
  return(X)
  
}, X = month_split)

rain_moon <- bind_rows(month_split) %>%
  tibble() %>%
  select(site_id, visit, precipitation, moon_fraction)

# Detection covariates combined -------------------------------------------

detection_covariates <- left_join(bind_rows(sites_in_grid), rain_moon, by = c("unique_site" = "site_id", "visit")) %>%
  left_join(trap_nights, by = c("unique_site", "visit")) %>%
  arrange(unique_site) %>%
  rename(site_id = unique_site)

write_rds(detection_covariates, here("data", "processed_data", "detection_covariates.rds"))


# Occurrence covariates ---------------------------------------------------
# The primary outcome is the effect of habitat type on occurrence we extract this from the site data

if(!file.exists(here("data", "processed_data", "occurrence_covariates.rds"))) {
  
  land_use <- bind_rows(sites_in_grid) %>%
    select(site_id = unique_site, landuse)
  
  # Get the bounding box of the village and buffer it by 100m before downloading from OSM
  
  distance_from_building <- function(data = spatial_traps, village_name) {
    
    osm <- opq(st_as_sfc(st_bbox(data %>%
                                   filter(village == village_name) %>%
                                   st_as_sf(coords = c("site_easting", "site_northing"), crs = SL_UTM) %>%
                                   st_transform(crs = default_CRS))) %>%
                 st_buffer(dist = 100) %>%
                 st_bbox()) %>%
      add_osm_feature(key = "building") %>%
      osmdata_sf()
    
    buildings <- osm$osm_polygons %>%
      st_transform(crs = SL_UTM) %>%
      st_union()
    
    sites <- bind_rows(sites_in_grid) %>%
      filter(village == village_name) %>%
      distinct(site_id = unique_site, site_easting, site_northing) %>%
      st_as_sf(coords = c("site_easting", "site_northing"), crs = SL_UTM) %>%
      mutate(distance_building = as.numeric(st_distance(., buildings)))
  }
  
  distance_building <- lapply(X = c("baiama", "lalehun", "lambayama", "seilama"), function(X) {
    distance_from_building(village_name = X)
  }) %>%
    bind_rows() %>%
    tibble() %>%
    distinct(site_id, distance_building)
  
  # We also use the distance from the centre of the village site these are stored as coordinates
  village_coords <- tibble(village = c("baiama", "lalehun", "lambayama", "seilama"),
                           X = c(-11.268454, -11.0803, -11.198249, -11.193628469657279),
                           Y = c(7.83708, 8.197533, 7.854131, 8.122285428353395)) %>%
    st_as_sf(coords = c("X", "Y"), crs = default_CRS) %>%
    st_transform(crs = SL_UTM)
  
  distance_from_centre <- bind_rows(sites_in_grid) %>%
    distinct(site_id = unique_site, site_easting, site_northing) %>%
    st_as_sf(coords = c("site_easting", "site_northing"), crs = SL_UTM) %>%
    mutate(distance_centre = case_when(str_detect(site_id, "baiama") ~ as.numeric(st_distance(., village_coords %>%
                                                                                                filter(village == "baiama"))),
                                       str_detect(site_id, "lalehun") ~ as.numeric(st_distance(., village_coords %>%
                                                                                                 filter(village == "lalehun"))),
                                       str_detect(site_id, "lambayama") ~ as.numeric(st_distance(., village_coords %>%
                                                                                                   filter(village == "lambayama"))),
                                       str_detect(site_id, "seilama") ~ as.numeric(st_distance(., village_coords %>%
                                                                                                 filter(village == "seilama"))),
                                       TRUE ~ as.numeric(NA))) %>%
    tibble() %>%
    distinct(site_id, distance_centre)
  
  # Elevation will also be used as an occurrence covariate
  # This method is currently failing due to an invalid/expired certificate
  # elevation_3s(lon = tile_coords[1], lat = tile_coords[2], path = here("data", "geodata"))
  # Will use the elevatr package instead for this we need a data.frame with the centre of each village
  
  elevation_rast <- get_elev_raster(locations = village_coords %>%
                                      st_transform(crs = default_CRS), prj = default_CRS, z = 12) %>%
    rast()
  
  elevation_vect <- vect(bind_rows(sites_in_grid) %>%
                           distinct(site_id = unique_site, site_easting, site_northing) %>%
                           st_as_sf(coords = c("site_easting", "site_northing"), crs = SL_UTM) %>%
                           st_transform(crs = default_CRS))
  
  elevation_vect$elevation <- terra::extract(elevation_rast, elevation_vect)[, 2]
  
  elevation <- data.frame(elevation_vect)
  
  population_rast <- population(2020, res = 0.5, path = here("data", "geodata"))
  
  population_vect <- elevation_vect[, 1]
  
  population_vect$population <- terra::extract(population_rast, population_vect)[, 2]
  
  population_vect$pop_quartile <- factor(ntile(population_vect$population, 4), levels = c(1, 2, 3, 4), labels = c("first", "second", "third", "fourth"))
  
  # Effectively this splits villages into different quartiles so this can be manually coded instead of the overlap that happens due to non-even split
  as.data.frame(population_vect) %>%
    mutate(village = str_split(site_id, "_", simplify = TRUE)[, 1]) %>%
    janitor::tabyl(village, pop_quartile)
  
  population_vect$pop_quartile <- as.data.frame(population_vect) %>%
    mutate(quartile = factor(case_when(str_detect(site_id, "seilama") ~ 1,
                                       str_detect(site_id, "baiama") ~ 2,
                                       str_detect(site_id, "lalehun") ~ 3,
                                       str_detect(site_id, "lambayama") ~ 4),
                             levels = c(1, 2, 3, 4), labels = c("first", "second", "third", "fourth"))) %>%
    pull(quartile)
  
  population <- data.frame(population_vect)
  
  # Occurrence covariates combined ------------------------------------------
  
  occurrence_covariates <- bind_rows(sites_in_grid) %>%
    distinct(site_id = unique_site, village) %>%
    left_join(land_use, by = c("site_id")) %>%
    left_join(distance_building, by = c("site_id")) %>%
    left_join(distance_from_centre, by = c("site_id")) %>%
    left_join(elevation, by = c("site_id")) %>%
    left_join(population, by = c("site_id")) %>%
    distinct() %>%
    arrange(site_id)
  
  write_rds(occurrence_covariates, here("data", "processed_data", "occurrence_covariates.rds"))
  
} else {
  
  occurrence_covariates <- read_rds(here("data", "processed_data", "occurrence_covariates.rds"))
  
}

# Site coordinates --------------------------------------------------------

coords <- bind_rows(sites_in_grid) %>%
  distinct(site_id = unique_site, site_easting, site_northing) %>%
  arrange(site_id)

coord_array <- array(data = c(coords$site_easting, coords$site_northing), dim = c(nrow(coords), 2), dimnames = list(c(coords$site_id), c("X", "Y")))

write_rds(coord_array, here("data", "processed_data", "site_coords.rds"))
