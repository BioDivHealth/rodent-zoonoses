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

# find enclosing coordinate boxes for files

#df <- host_path %>%
 # rowwise() %>%
#  mutate(Upper_Left_Corner = list(get_dec_box(decimalLongitude, decimalLatitude))) %>%
#  unnest_wider(Upper_Left_Corner, names_sep = "_") %>%
#  rename(Upper_Left_Lat = Upper_Left_Corner_1, Upper_Left_Lon = Upper_Left_Corner_2) %>%
#  dplyr::select(decimalLatitude, decimalLongitude, Upper_Left_Lat, Upper_Left_Lon)

# find unique boxes

#files <- unique(df[,c('Upper_Left_Lat','Upper_Left_Lon')])

# 3. Plot example raster from our GLC data
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

# map values to colours using legend for correct land cover colours
# Load the CSV file
legend_file <- "./data/example_raster/LC_legend.csv"
legend <- read.csv(legend_file)

legend$colour <- sapply(legend$colour, parse_rgb)

# Create a lookup table from the legend file
lookup_table <- legend %>%
  dplyr::select(number, colour)

# Convert the raster to a data frame for ggplot2
r_df <- as.data.frame(cropped_ras, xy = TRUE, na.rm = TRUE)

# Map the raster values to colors
r_df <- r_df %>%
  left_join(lookup_table, by = c("X19852000_E15N45_3" = "number"))

# add buffer to plot
# Define the center point for the buffer circle (longitude, latitude)
center_point <- st_point(c(16.56, 44.06))

# Convert to an sf object
center_sf <- st_sfc(center_point, crs="+proj=longlat +datum=WGS84 +no_defs")

# Define the buffer radius 
buffer_radius <- 1000  # 5 km radius

# Create the buffer for ggplot
buffer_circle <- st_buffer(center_sf, dist = buffer_radius)

# Convert the buffer circle to a data frame for ggplot2
buffer_df <- as.data.frame(st_coordinates(st_cast(buffer_circle, "POLYGON")))

# plot on raster with ggplot

ggplot(r_df, aes(x = x, y = y, fill = colour)) +
  geom_raster() +
  scale_fill_identity() +  # Use the exact colors provided
  geom_polygon(data = buffer_df, aes(x = X, y = Y), fill = "red", alpha=0.35) +
  theme_minimal() +
  labs(fill = "Colour", x = "Longitude", y = "Latitude") +
  geom_point(x=16.56, y=44.06, color = "red")
  coord_fixed()

ggsave('./figures/bufferplot.pdf',device="pdf", width = 9, height = 9)

# Extract the values within the buffer
land_cover_values <- raster::extract(cropped_ras, buffer)

# Flatten the list of values
land_cover_values <- unlist(land_cover_values)

# Convert numbers to habitat names using legend csv
# Load the CSV file
dict_file <- "./data/example_raster/LC_legend.csv"
dict <- read.csv(dict_file)

# Create a lookup table
lookup_table2 <- dict %>%
  dplyr::select(number, habitat)

lookup_table3 <- legend %>%
  dplyr::select(habitat, colour)

lookup_vector <- setNames(lookup_table2$habitat, lookup_table2$number)

# Convert the numbers
converted_values <- lookup_vector[as.character(land_cover_values)]


# Calculate the percentage of each land cover type
land_cover_table <- table(converted_values)
land_cover_percentage <- 100 * land_cover_table / sum(land_cover_table)
land_cover_df <- data.frame(land_cover_percentage)

land_cover_df <- land_cover_df %>%
  left_join(lookup_table3, by = c("converted_values" = "habitat"))

#make dictionary for correct legend names
legend_dict <- setNames(as.character(lookup_table3$habitat), as.character(lookup_table3$colour))

# Plot the results as pie chart
ggplot(land_cover_df, aes(x = "", y = Freq, fill=colour)) +
  scale_fill_identity(guide="legend", labels= legend_dict) +
  theme_minimal() +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  labs(title = "Percentage of land cover", fill = "Land cover type")
