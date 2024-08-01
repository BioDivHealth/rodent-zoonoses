# 1. Load packages ----
pacman::p_load(here,
               stringr,
               dplyr,
               ggplot2,
               raster,
               tiff)  

source(here("scripts", "00_useful_functions.R")) #loaded as list `myfuncs`

# 2. Load example GLC data file
  
str_name<-'./data/GLC/GLC_FCS30D_20002022_W175N70_Annual.tif' 
eg_raster=raster(str_name)

eg_raster <- setMinMax(eg_raster)

plot(eg_raster, main="example of GLC_FCS raster")

# 3. load our data to find necessary files to download

host_path <- readRDS("data/host_path_wide.rds")

files <- glc.files(host_path, "decimalLongitude")
files




