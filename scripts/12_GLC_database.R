# 1. Load packages ----
pacman::p_load(here,
               stringr,
               dplyr,
               ggplot2,
               raster,
               tiff,
               tidyr,
               sp,
               rasterVis)  

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

### make raster stack
str_name <-'./data/example_raster/19852000_E15N45.tif' 
raster_stack <- stack(str_name)

# select 1990 as example
croat_1990 <- raster_stack[[2]]

plot(croat_1990, main="Land cover for Croatia in 1990")

# 4. find percentage land cover around coordinates

### find example coordinate in our data

# Define the ranges for filtering
long_min <- 15
long_max <- 20
lat_min <- 40
lat_max <- 45

# Filter the dataset
filtered_df <- host_path %>%
  filter(decimalLongitude >= long_min & decimalLongitude <= long_max,
         decimalLatitude >= lat_min & decimalLatitude <= lat_max)

latitude <- filtered_df$decimalLatitude[1]
longitude <- filtered_df$decimalLongitude[1]
year <- floor_to_bin(as.numeric(format(filtered_df$start_date[1], "%Y")))

# Create dict to convert year to slice of raster stack
dict <- c("1985" = 1, "1990" = 2, "1995" = 3)

# Convert the number to a character string to match the named vector keys
slice <- dict[as.character(year)]

croat_1995 <- raster_stack[[slice]]

# Define the buffer radius
buffer_radius <- 1000  # in meters

# Create a spatial point for the coordinate
coord <- SpatialPoints(cbind(longitude, latitude), 
                       proj4string = CRS(projection(croat_1995)))

# Create a buffer around the coordinate
buffer <- buffer(coord, width = buffer_radius)

# Plot the buffer and the land cover raster for visualization

#crop raster
b <- as(extent(16.53, 16.59, 44.03, 44.07), 'SpatialPolygons')
crs(b) <- crs(croat_1995)
cropped_ras <- crop(croat_1995, b)

plot(cropped_ras)
plot(buffer, add = TRUE)

# coord_cartesian() within ggplot for zooming in plot area

# Extract the values within the buffer
land_cover_values <- extract(croat_1995, buffer)

# Flatten the list of values
land_cover_values <- unlist(land_cover_values)

# Calculate the percentage of each land cover type
land_cover_table <- table(land_cover_values)
land_cover_percentage <- 100 * land_cover_table / sum(land_cover_table)

# Print the results
print(land_cover_percentage)




















