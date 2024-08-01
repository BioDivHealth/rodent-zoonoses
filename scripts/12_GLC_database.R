# 1. Load packages ----
pacman::p_load(here,
               stringr,
               dplyr,
               ggplot2,
               raster,
               tiff,
               tidyr)  

source(here("scripts", "00_useful_functions.R")) #loaded as list `myfuncs`

# 2. load our data to find necessary files to download

host_path <- readRDS("data/host_path_wide.rds")

# find enclosing coordinate boxes for files

df <- host_path %>%
  rowwise() %>%
  mutate(Upper_Left_Corner = list(get_dec_box(decimalLongitude, decimalLatitude))) %>%
  unnest_wider(Upper_Left_Corner, names_sep = "_") %>%
  rename(Upper_Left_Lat = Upper_Left_Corner_1, Upper_Left_Lon = Upper_Left_Corner_2) %>%
  dplyr::select(decimalLatitude, decimalLongitude, Upper_Left_Lat, Upper_Left_Lon)

# find unique boxes

files <- unique(df[,c('Upper_Left_Lat','Upper_Left_Lon')])

# 3. Plot example raster from our GLC data
str_name <-'./data/example_raster/19852000_E15N45.tif' 
eg_raster=raster(str_name)

eg_raster <- setMinMax(eg_raster)

plot(eg_raster, main="example of GLC_FCS raster")
