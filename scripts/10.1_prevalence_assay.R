library(dplyr)
library(here)
library(stringr)
library(ggplot2)
library(patchwork)
library(gridExtra)
library(viridis)

pathogen1 <- readRDS("data/pathogen_assay_record.rds")

pathogen <- pathogen1 %>% 
  mutate(
    pathogen_abbrev = case_when(
      scientificName == "Sin Nombre Virus" ~ "SNV",
      scientificName == "Sin Nombre Orthohantavirus" ~ "SNV", 
      scientificName == "Sin Nombre virus" ~ "SNV", 
      scientificName == "Catarina Virus" ~ "CTNV",
      scientificName == "Oliveros Virus" ~ "OLV",
      scientificName == "Latino mammarenavirus" ~ "LATV",
      scientificName == "Orthohantavirus dobravaense" ~ "DOBV",
      scientificName == "Dobrava-Belgrade orthohantavirus" ~ "DOBV",
      scientificName == "Mammarenavirus bearense" ~ "BCNV",
      scientificName == "Puumala orthohantavirus" ~ "PUUV",
      scientificName == "Puumala virus" ~ "PUUV",
      scientificName == "Seoul virus" ~ "SEOV",
      scientificName == "Guanarito virus" ~ "GTOV",
      scientificName == "Lassa mammarenavirus" ~ "LASV",
      scientificName == "Pirital virus" ~ "PIRV",
      scientificName == "Rat hepatitits E virus" ~ "RHEV",
      scientificName == "Hantaan virus" ~ "HTNV",
      scientificName == "Hantaviridae" ~ "HTNV",
      scientificName == "mammarenavirus" ~ "",
      scientificName == "Leptospriosa interrogans" ~ "",
      # Keep all other names the same
      TRUE ~ scientificName
    ),
  )

# Gives rodents tested >1 time for same pathogen using different assays
rodent_test_counts <- pathogen %>%
  # Convert numeric variables to numeric class
  mutate(
    positive = as.numeric(positive),
    tested = as.numeric(tested),
    
  ) %>%
  # Get rid of rows where 0 tested, where tested is NA, where positives are NA
  filter(tested != 0, !is.na(tested), !is.na(positive)) %>%
  group_by(associated_rodent_record_id, pathogen_abbrev) %>%
  filter(n_distinct(assay) > 1) %>%
  ungroup() %>%
  mutate(prevalence = positive / tested) %>%
  arrange(associated_rodent_record_id)

# Find unique assay names
rodent_test_counts$assay %>% unique()

# simplify naming of assays by categorising into broad groups
rodent_assay_counts <- rodent_test_counts %>%
  mutate(simple.assay = 
           case_when(
             str_detect(assay, "ELISA") | str_detect(assay, "EL") | str_detect(assay, "erology") | str_detect(assay, "Immuno")| str_detect(assay, "IFA") | str_detect(assay, "IB") ~ "Serology",
             str_detect(assay, "PCR")   ~ "PCR",
           ))

ggplot(rodent_assay_counts, aes(x = prevalence, y = associated_rodent_record_id, group = simple.assay)) +
  geom_point(aes(color = simple.assay, size = tested)) +  # Dots representing each pathogen
  geom_line(aes(group = associated_rodent_record_id), color = "grey", size = 1) +  # Lines connecting dots
  labs(
    title = "Dumbbell Plot of Rodents with Multiple Pathogen Tests",
    x = "Prevalence",
    y = "Rodent Record ID",
    color = "Assay"
  ) +
  theme_minimal()


# Fixing graph
# Create different plots for each pathogen with rodent_id with >1 tests
all_pathogens_plot <- ggplot(rodent_assay_counts, aes(x = prevalence, y = associated_rodent_record_id, group = simple.assay)) +
  geom_point(aes(color = simple.assay, size = tested), alpha = 0.7) +  # Adjust size of points
  geom_line(aes(group = associated_rodent_record_id), color = "grey", size = 0.5) +  # Thinner lines
  labs(
    title = "Rodent_IDs with Multiple Pathogen Tests",
    x = "Prevalence",
    y = "Rodent Record ID",
    color = "Assay"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 8),  # Smaller text
    axis.title = element_text(size = 9),  # Axis title size
    axis.text = element_text(size = 7),  # Axis text size
    legend.title = element_text(size = 8),  # Legend title size
    legend.text = element_text(size = 7)  # Legend text size
  ) +
  facet_wrap(~ scientificName, scales = "free_y")

# Create a separate plot specifically for Sin Nombre Orthohantavirus with smaller text and points
sin_nombre_plot <- ggplot(rodent_assay_counts %>% filter(scientificName == "Sin Nombre Orthohantavirus"), aes(x = prevalence, y = associated_rodent_record_id, group = simple.assay)) +
  geom_point(aes(color = simple.assay, size = tested), alpha = 0.7) +  # Adjust size of points
  geom_line(aes(group = associated_rodent_record_id), color = "grey", size = 0.5) +  # Thinner lines
  labs(
    title = "Sin Nombre Orthohantavirus",
    x = "Prevalence",
    y = "Rodent Record ID",
    color = "Assay"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 8),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7)
  )

# Visualize them
print(all_pathogens_plot)
print(sin_nombre_plot)

# Combine them vertically
grid.arrange(all_pathogens_plot, sin_nombre_plot, ncol = 1)



### Calculate density vs Prevalence

# PCR plot
pcr_plot <- ggplot(subset(rodent_assay_counts, simple.assay == "PCR"), 
                   aes(x = prevalence, fill = pathogen_abbrev)) +
  geom_density(alpha = 0.5, bw = 0.1) +  # Density plot
  labs(
    title = "Density of Prevalence for PCR Assays",
    x = "Prevalence",
    y = "Density",
    fill = "Pathogen"  # Legend title
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 8),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7)
  )

# Serology plot
serology_plot <- ggplot(subset(rodent_assay_counts, simple.assay == "Serology"), 
                        aes(x = prevalence, fill = pathogen_abbrev)) +
  geom_density(alpha = 0.5, bw = 0.1) +  # Density plot
  labs(
    title = "Density of Prevalence for Serology Assays",
    x = "Prevalence",
    y = "Density",
    fill = "Pathogen"  # Legend title
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 8),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7)
  )

# Display the plots together
library(gridExtra)
grid.arrange(pcr_plot, serology_plot, ncol = 1)

#Bar plot without grouping by rodent_ID
ggplot(rodent_assay_counts, aes(x = pathogen_abbrev, y = prevalence, fill = simple.assay)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = c("PCR" = "darkblue", "Serology" = "darkgreen")) +
  labs(title = "Prevalence by Assay Type",
       x = "Pathogen",
       y = "Prevalence",
       fill = "Assay") +
  theme_bw()

#Density graph without grouping by rodent_ID
ggplot(rodent_assay_counts, aes(x = prevalence, fill = pathogen_abbrev)) +
  geom_density(alpha = 0.6, bw = 0.05) +  # Density plot with transparency
  facet_wrap(~ simple.assay, scales = "free") +  # Separate plots for PCR and Serology
  labs(title = "Density of Prevalence by Pathogen and Assay Type",
       x = "Prevalence",
       y = "Density",
       fill = "Pathogen") +
  theme_bw() +
  scale_fill_viridis_d(option = "plasma")  # Use the 'plasma' option, or try 'viridis', 'magma', etc.



