library(readr)
library(dplyr)
install.packages("writexl")
library(writexl)
install.packages("readxl")
library(readxl)

# Create list from Dave cleaned data (data_for_ana) of all host and pathogens conserving uniq_ids
# Create aggregated list of host and pathogens without uniq_ids and summing ttested and positve for rows with same virus and species
# Read the CSV file (assuming your file is in the current working directory)
file_path <- "~/Desktop/pathogen_host list.xlsx"

# Read the Excel file, specifying the sheet name if needed
data_uniq_id <- read_excel(file_path, sheet= "Sheet1")

# Filter out rows where the 'tested' column has the value 'NA'
cleaned_data_uniq_id <- data_uniq_id %>% filter(tested != "NA")

# Filter out rows where the 'positive' column has the value '0'
cleaned_data2_uniq_id <- data_uniq_id %>%
  filter(tested != "NA" & positive != 0)

write_xlsx(cleaned_data2_uniq_id, "~/Desktop/host-pathogenlist2.xlsx")



# Create aggregated list of host and pathogens without uniq_ids and sum up tested and positive for rows with same virus and species
# Read the CSV file
file_path <- "~/Desktop/pathogen_host list.xlsx"

# Read the Excel file, specifying the sheet name if needed
data <- read_excel(file_path, sheet= "Sheet2")

# Filter out rows where the 'tested' column has the value 'NA'
cleaned_data <- data %>% filter(tested != "NA")

# Filter out rows where the 'positive' column has the value '0'
cleaned_data2 <- data %>%
  filter(tested != "NA" & positive != 0)


# Aggregate repetitive rows (considering only word variables) and sum up number variables from repetitive rows
# List of word columns
word_columns <- c("family", "order", "pathogen_genus...3", "pathogen_family...4", "pathogen_taxa...5", "pathogen_scientificName", "pathogen_genus...9", "pathogen_family...10", "pathogen_taxa...11")
# List of numeric columns
numeric_columns <- c("tested", "positive")

# Clean the word columns to remove any extra spaces and ensure consistent formatting
data2 <- cleaned_data2 %>%
  mutate(across(all_of(word_columns), ~trimws(as.character(.))))

# Ensure the numeric columns are numeric
data2 <- cleaned_data2 %>%
  mutate(across(all_of(numeric_columns), as.numeric))

# Aggregate the data by summing numeric columns for each group of word columns
aggregated_data <- data2 %>%
  group_by(across(all_of(word_columns))) %>%
  summarise(across(all_of(numeric_columns), sum, .names = "{.col}")) %>%
  ungroup()

# Save data in desktop
write_xlsx(aggregated_data, "~/Desktop/host-pathogenlist.xlsx")
