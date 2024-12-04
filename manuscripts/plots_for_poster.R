########  PLOT 1

# Load necessary libraries
library(ggplot2)
library(viridis)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(tidyverse)

# Load data frame
host_path <- readRDS("../data/host_path_wide.rds")
attach(host_path)

# Get the world map data
world <- ne_countries(scale = "large", returnclass = "sf")
coords_sf <- st_as_sf(host_path, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326, agr = "constant")
host_path$year <- as.numeric(format(host_path$start_date, "%Y"))

# Plot the points onto the world map 
ggplot(data = world) +
  geom_sf(colour="darkgrey") +
  geom_sf(data = coords_sf, aes(geometry = geometry, colour=host_path$year), size = 3, shape=19) +
  coord_sf(expand = FALSE) +                    # Keep the aspect ratio
  theme_void() +
  theme(panel.background = element_rect(), 
        axis.text = element_blank(),            # Remove axis text (tick labels)
        axis.ticks = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "cm")) +
  scale_colour_viridis(name="Year", direction=1, option="inferno", limits=c(1980, 2030))

ggsave("map_plot.png", dpi=1200)

########  PLOT 2

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(scico)

# group by host and sum count
host_plot <- host_path %>% group_by(genus) %>% summarise(n_hosts=sum(as.numeric(individualCount)))
host_plot <- host_plot[-1, ] %>% arrange(desc(n_hosts)) %>% slice(1:12)

# Plot barchart
ggplot(host_plot, aes(x = reorder(genus, -n_hosts), y = n_hosts, fill=genus)) +
  geom_bar(stat = "identity") + 
  annotate("point", x = 9, y = 1500, shape = 21, size = 50, fill = "chocolate4", color = "chocolate4") +
  annotate("text", x = 9, y = 1500, label = "10955 rodents\n total", size = 5) +
  labs(x = "Genus",
       y = "Number Trapped") +
  theme_minimal() +
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),
        axis.text.x = element_text(angle = 45, hjust = 1, size=11),
        legend.title = element_blank(),
        legend.position = "none") +
  scale_fill_scico_d(palette = "tokyo")

ggsave("rodent_plot.png", dpi=1200)


########  PLOT 3

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(scico)

host_path <- readRDS("../data/host_path_wide.rds")
attach(host_path)

# group by virus and sum infections
infection_plot <- host_path %>% group_by(pathogen_abbrev) %>% summarise(n_infections=sum(positive))
infection_plot <- infection_plot[-1, ] %>% arrange(desc(n_infections))

# Plot barchart
ggplot(infection_plot, aes(x = reorder(pathogen_abbrev, -n_infections), y = n_infections, fill=pathogen_abbrev)) +
  geom_bar(stat = "identity") + 
  annotate("point", x = 10.45, y = 70, shape = 21, size = 50, fill = "cadetblue", color = "cadetblue") +
  annotate("text", x = 10.45, y = 70, label = "480 positive\n rodents", size = 5) +
  labs(x = "Pathogen",
       y = "Number of Infections") +
  theme_minimal() +
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),
        axis.text.x = element_text(size=11),
        legend.title = element_blank(),
        legend.position = "none") +
  scale_fill_scico_d(palette = "roma")

ggsave("pathogen_plot.png", dpi=1200)
