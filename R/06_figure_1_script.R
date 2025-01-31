source(here::here("R", "00_setup.R"))

# A map of Sierra Leone in Africa -----------------------------------------

world <- ne_countries(scale = "medium", returnclass = "sf")

africa <- world %>%
  filter(str_detect(continent, "Africa")) %>%
  select(admin) %>%
  mutate(fill = case_when(str_detect(admin, "Sierra Leone") ~ "black",
                          TRUE ~ "grey")) %>%
  ms_filter_islands(min_area = 1e10)

sl_bbox <- africa %>%
  filter(str_detect(admin, "Sierra Leone")) %>%
  st_bbox() %>%
  st_as_sfc()

africa_map <- ggplot() +
  geom_sf(data = africa, aes(fill = fill)) +
  geom_sf(data = sl_bbox, fill = NA, colour = "black", size = 4) +
  scale_fill_manual(values = c("grey", "white")) +
  guides(fill = "none") +
  theme_void()

sle_sf <- geodata::gadm(country = "SLE", level = 2, path = here("data", "geodata")) %>%
  st_as_sf()

fig_1_palette <- c(village_palette, "#FFFFFF")
names(fig_1_palette) <- c(names(village_palette)[1:4], "poi")

poi <- tibble(name = c("Kenema", "Baiama", "Lalehun", "Lambayama", "Seilama"),
              cat = c("poi", "Baiama", "Lalehun", "Lambayama", "Seilama"),
              lat = c(7.876161956810467, 7.837529372181356, 8.197392257077409, 7.850593096948891, 8.12230048178563),
              lon = c(-11.190811585001954, -11.268407665149846, -11.08032958100431, -11.196939025872055, -11.1935976318381)) %>%
  st_as_sf(coords = c("lon", "lat"), crs = default_CRS)


# Focus map on Kenema
kenema_map <- get_googlemap(st_coordinates(poi %>% filter(name == "Kenema"))[1, ], zoom = 9, maptype = "satellite", scale = 2) %>%
  rast() %>%
  mask(vect(sle_sf %>%
              filter(str_detect(NAME_2, "Kenema"))))

zoom_eastern <- sle_sf %>%
  filter(str_detect(NAME_1, "Eastern")) %>% st_bbox()

# sle_map <- ggplot() +
#   geom_sf(data = sle_sf) +
#   geom_sf(data = sle_sf %>%
#             filter(str_detect(NAME_2, "Kenema")), fill = "orange") +
#   geom_sf(data = zoom_kenema %>%
#             st_as_sfc(), fill = NA, colour = "black", size = 4)  +
#   guides(fill = "none") +
#   theme_void()

study_map <- ggplot() + 
  geom_sf(data = sle_sf, fill = "grey", colour = "white") +
  geom_sf(data = sle_sf %>%
            filter(str_detect(NAME_1, "Eastern")), fill = NA, colour = "black") +
  geom_sf(data = poi, size = 2, aes(shape = name)) +
  ggrepel::geom_text_repel(data = poi, aes(label = name, geometry = geometry), stat = "sf_coordinates", min.segment.length = 0) +
  labs(x = element_blank(),
       y = element_blank()) +
  theme_bw() +
  scale_shape_manual(values = c(17, 19, 18, 15, 3)) +
  coord_sf(xlim = c(zoom_eastern[1], zoom_eastern[3]), ylim = c(zoom_eastern[2], zoom_eastern[4])) +
  ggspatial::annotation_north_arrow(style = north_arrow_minimal()) +
  ggspatial::annotation_scale(location = "br") +
  guides(shape = "none") +
  labs(title = "A)")

sl_inset_map <- ggdraw() +
  draw_plot(study_map) +
  draw_plot(africa_map, x = 0.7, y = 0.6, width = 0.3, height = 0.3)

# Raster version ----------------------------------------------------------
library("tidyterra")
SLE_rast <- rast(here("data", "processed_data", "SLE_habitat.tif"))
site_coords <- read_rds(here("data", "processed_data", "site_coords.rds"))

site_coords_vect <- site_coords %>%
  vect(crs = SL_UTM)
site_coords_vect$site = str_to_sentence(str_split(row.names(site_coords), "_", simplify = TRUE)[, 1])
site_poly <- site_coords_vect %>%
  convHull(by = "site")

iucn_clr <- read_csv(here("data", "geodata", "iucn_clr.csv"), col_names = FALSE) %>%
  data.frame()
coltab(SLE_rast) <- iucn_clr[1:5]
levels(SLE_rast) <- iucn_clr[c(1,6)]

sle_sf_2 <- geodata::gadm(country = "SLE", level = 2, path = here("data", "geodata")) %>%
  st_as_sf()

eastern_vect <- sle_sf_2 %>%
  filter(str_detect(NAME_2, "Kenema")) %>%
  vect()

crop_rast <- crop(SLE_rast, eastern_vect)

sle_map_2 <- ggplot() +
  geom_sf(data = sle_sf_2, fill = "grey") +
  geom_sf(data = sle_sf_2 %>%
            filter(NAME_2 == "Kenema"), fill = "orange") +
  theme_nothing()

study_map_2 <- ggplot() +
  geom_spatraster(data = crop_rast, use_coltab = FALSE) +
  geom_spatvector(data = eastern_vect, fill = NA, colour = "black", lwd = 1) +
  geom_sf(data = poi, size = 2, fill = "black", aes(shape = name)) +
  scale_shape_manual(values = c(17, 19, 18, 15, 3)) +
  ggrepel::geom_label_repel(data = poi, aes(label = name, geometry = geometry), stat = "sf_coordinates", min.segment.length = 0, size = 4) +
  scale_fill_wiki_d() +
  labs(fill = "IUCN Habitat",
       x = element_blank(),
       y = element_blank(),
       shape = element_blank()) +
  guides(shape = "none") +
  ggspatial::annotation_north_arrow(style = north_arrow_minimal()) +
  ggspatial::annotation_scale(location = "br") +
  labs(title = "A)") +
  theme_bw() + 
  theme(legend.text = element_text(size = 12))

sl_inset_map_v2 <- ggdraw() +
  draw_plot(study_map_2) +
  draw_plot(sle_map_2, x = 0.7, y = 0.7, width = 0.3, height = 0.3)

ggsave2(plot = sl_inset_map_v2, filename = here("output", "figures", "Figure_1a.png"), dpi = 300, width = 8, height = 6)

# Trap timeline -----------------------------------------------------------
trap_data <- read_csv(here("data", "input", "trap_data.csv"))

timeline <- trap_data %>%
  tibble() %>%
  select(-geometry) %>%
  filter(village != "bambawo") %>%
  group_by(village, visit, grid_number, trap_number) %>%
  summarise(tn = n(),
            date_set = min(date_set)) %>%
  group_by(visit, village) %>%
  summarise(tn = sum(tn),
            date_set = min(date_set)) %>%
  ungroup() %>%
  mutate(village = str_to_title(village),
         date_set = case_when(date_set == as.Date("2023-02-05") ~ as.Date("2023-02-08"),
                              TRUE ~ date_set))

rain_season <- tibble(season = "rainy",
                      date_start = c(as_date("2021-05-01"), as_date("2022-05-01"), as_date("2023-05-01")),
                      date_end = c(as_date("2021-11-01"), as_date("2022-11-01"), as_date("2023-11-01")))

timeline_plot <- ggplot(timeline) +
  geom_rect(data = rain_season, aes(xmin = date_start, xmax = date_end, ymin = 0, ymax = Inf), fill = "lightblue", alpha = 0.5) +
  geom_point(aes(x = date_set, y = tn, shape = village), size = 1) +
  coord_cartesian(xlim = c(as.Date("2020-11-01"), as.Date("2023-05-01"))) +
  scale_shape_manual(values = c(17, 18, 15, 3)) +
  theme_bw() +
  labs(x = "Visit date",
       y = "Trap nights",
       shape = "Village") +
  labs(title = "B)") +
  theme(legend.text = element_text(size = 12))

save_plot(plot = timeline_plot, filename = here("output", "figures", "Figure_1b.png"), base_height = 3, base_width = 8)

save_plot(plot = plot_grid(plotlist = list(sl_inset_map_v2,
                                           timeline_plot),
                           ncol = 1,
                           rel_heights = c(8, 2)),
          filename = here("output", "figures", "Figure_1_combined.png"),
          base_height = 10,
          base_width = 8)
