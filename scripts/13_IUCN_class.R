### This script loads the IUCN land cover classification legend. We hope to use some method e.g. decision tree to 
### identify equivalent land cover classifications in the GLC database.

# 1. Load packages ----
pacman::p_load(here,
               stringr,
               dplyr,
               ggplot2,
               raster,
               tiff,
               tidyr,
               sf)  

source(here("scripts", "00_useful_functions.R")) #loaded as list `myfuncs`

# 2. load our data to find necessary files to download

host_path <- readRDS("data/host_path_wide.rds")

str_name <-'./data/IUCN/iucn_habitatclassification_composite_lvl2_ver003.tif' 
raster_stack <- stack(str_name)


plot(raster_stack, main="IUCN LC class")
