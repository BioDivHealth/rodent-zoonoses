# 1. Load packages ----
pacman::p_load(here,
               stringr,
               dplyr,
               ggplot2,
               raster,
               tiff)  

# 2. Load example GLC data file



str_name<-'./GLC/GLC_FCS30D_20002022_W175N70_Annual.tif' 
eg_raster=raster(str_name)

eg_raster <- setMinMax(eg_raster)

plot(eg_raster, main="example of GLC_FCS raster")
