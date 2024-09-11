library(dplyr)
library(here)
library(stringr)

pathogen <- readRDS("data/pathogen_assay_record.rds")

studies <- subset(combined_data, studies == "site-session")

# Gives rodents tested >1 time for same pathogen using different assays
rodent_test_counts <- pathogen %>%
  # Convert numeric variables to numeric class
  mutate(
    positive = as.numeric(positive),
    tested = as.numeric(tested),
    
  ) %>%
  # Get rid of rows where 0 tested, where tested is NA, where positives are NA
  filter(tested != 0, !is.na(tested), !is.na(positive)) %>%
  group_by(associated_rodent_record_id, scientificName) %>%
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
             str_detect(assay, "ELISA") | str_detect(assay, "erology") | str_detect(assay, "Immuno") ~ "Serology",
             str_detect(assay, "PCR")   ~ "PCR",
           ))

ggplot(rodent_assay_counts, aes(x = associated_rodent_record_id, y = prevalence, group = simple.assay)) +
  geom_point(aes(color = simple.assay, size = tested)) +  # Dots representing each pathogen
  geom_line(aes(group = associated_rodent_record_id), color = "grey", size = 1) +  # Lines connecting dots
  labs(
    title = "Dumbbell Plot of Rodents with Multiple Pathogen Tests",
    x = "Rodent Record ID",
    y = "Prevalence",
    color = "Assay"
  ) +
  theme_minimal()
