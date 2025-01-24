
# Install and load required packages --------------------------------------

if (!require("pacman")) install.packages("pacman")

pkgs =
  c("coda",
    "cowplot",
    "DescTools",
    "elevatr",
    "fastDummies",
    "flextable",
    "geodata",
    "ggmap",
    "ggnewscale",
    "ggridges",
    "ggspatial",
    "ggtext",
    "googledrive",
    "here",
    "lubridate",
    "mapview",
    "osmdata",
    "RColorBrewer",
    "RhpcBLASctl",
    "rmapshaper",
    "rnaturalearth",
    "rnaturalearthdata",
    "rosm",
    "sf",
    "spOccupancy",
    "stars",
    "suncalc",
    "terra",
    "tidyfast",
    "tidygraph",
    "tidyterra",
    "tidyverse",
    "vegan"
  )

pacman::p_load(pkgs, character.only = T)


# Define Coordinate Reference System for project --------------------------

default_CRS <- "EPSG:4326"
SL_UTM <- "EPSG:32629"


# Define order of villages for consistency --------------------------------

village_order <- c("baiama","lalehun", "lambayama", "seilama")

village_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#888888")
names(village_palette) <-  c("Lalehun", "Seilama", "Lambayama", "Baiama")

landuse_palette <- c("#00913A", "#FEC44F", "#A13B9E")
names(landuse_palette) <- c("Forest", "Agriculture", "Village")


# Define order of species for consistency ---------------------------------

# Species are ordered by n detected except the Shrew species is moved to after the Rodent species

all_species_order <- c("Mastomys natalensis", "Praomys rostratus", "Rattus rattus", "Mus musculus", "Lophuromys sikapusi", "Mus setulosus",
                       "Crocidura olivieri", "Crocidura buettikoferi", "Crocidura grandiceps", "Malacomys edwardsi", "Lemniscomys striatus",
                       "Hylomyscus simus", "Hybomys planifrons", "Mastomys erythroleucus", "Crocidura theresae", "Gerbilliscus guineae", "Dasymys rufulus")

species_order_plots <- c("Mastomys natalensis", "Praomys rostratus", "Rattus rattus", "Mus musculus", "Lophuromys sikapusi", "Mus setulosus",
                         "Crocidura olivieri", "Crocidura buettikoferi", "Crocidura grandiceps", "Malacomys edwardsi", "Lemniscomys striatus",
                         "Hylomyscus simus", "Hybomys planifrons", "Mastomys erythroleucus", "Crocidura theresae", "Gerbilliscus guineae", "Dasymys rufulus")


# Consistent colours for land use -----------------------------------------

group_landuse_palette <- c("#00913a", "#FEC44F", "#F7A820", "#A13B9E", "#5407A6")
names(group_landuse_palette) <- c("Forest - Rural", "Agriculture - Rural", "Agriculture - Peri-urban", "Village - Rural", "Village - Peri-urban")


# Define seasons associated with trapping sessions ------------------------

season <- tibble(visit = 1:10, season = c(rep("Dry", 2), rep("Rainy", 2), rep("Dry", 2), rep("Rainy", 2), rep("Dry", 2)))
