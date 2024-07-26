setwd("../data")

## Load packages
if (!require("pacman")) install.packages("pacman")

pkgs <- c(
  "here",
  "tidyverse",
  "lubridate",
  "dplyr"
)


pacman::p_load(pkgs, character.only = T)

## if multiple in put sheets, load and merge sheets
##harry
studies_h = read.csv("./aren_hant_data/studies_h.csv")
host_h = read.csv("./aren_hant_data/host_h.csv")
pathogen_h = read.csv("./aren_hant_data/pathogen_h.csv")

##david
studies_d = read.csv("./aren_hant_data/studies_d.csv")
host_d = read.csv("./aren_hant_data/host_d.csv")
pathogen_d = read.csv("./aren_hant_data/pathogen_d.csv")
pathogen_d = select(pathogen_d, -X, -X.1, -X.2, -X.3)

##ana
#studies_a = read.csv("./aren_hant_data/studies_a.csv")
#host_a = read.csv("./aren_hant_data/host_a.csv")
#host_a = select(host_a, -trapping.notes)
#pathogen_a = read.csv("./aren_hant_data/pathogen_a.csv")


studies <-  rbind(studies_h, studies_d)
host <-  rbind(host_h, host_d)
pathogen <-  rbind(pathogen_h, pathogen_d)

combined_data <- list(studies = studies,
                         host = host,
                         pathogen = pathogen)

write_rds(combined_data, file="combined_data_all.rds")
attach(combined_data)

## Filter for site resolution

studies <- subset(studies, data_resolution == "site-session")

# keep only host data with matching study ID from this subset

host <- host[host$study_id %in% studies$study_id, ]

## translate date into start/end columns

host <- host %>%
  separate(eventDate, into = c("start_date", "end_date"), sep = "/")

host_monthly <- host %>%
  filter(str_detect(start_date, "^\\d{4}-\\d{2}") | str_detect(start_date, "^\\d{4}-\\d{2}-\\d{2}"))

# Add "-01" to strings in "yyyy-mm" format
host_monthly$start_date <- ifelse(grepl("^\\d{4}-\\d{2}$", host_monthly$start_date), 
                            paste0(host_monthly$start_date, "-01"), 
                            host_monthly$start_date)

host_monthly$end_date <- ifelse(grepl("^\\d{4}-\\d{2}$", host_monthly$end_date), 
                                  paste0(host_monthly$end_date, "-01"), 
                                  host_monthly$end_date)

## convert all dates to yyyy-mm-dd and convert to date object

host_monthly$start_date <- as.Date(host_monthly$start_date)
host_monthly$end_date <- as.Date(host_monthly$end_date)

## Calculate duration of trapping
host_monthly$period <- host_monthly$end_date - host_monthly$start_date
host_monthly$period[is.na(host_monthly$period)] <- 0

## filter for less than 1 month trapping time
host <- host_monthly %>% filter(period <= 31)

## filter for multiple species
host <- host %>%
  group_by(study_id) %>%
  filter(n_distinct(scientificName) > 1) %>%
  ungroup()

## Filter pathogen dataset for presence of rodent_id in host dataset

pathogen <- pathogen[pathogen$associated_rodent_record_id %in% host$rodent_record_id, ]

## check for multiple pathogen entries for each rodent_id and take max positives

pathogen <- pathogen %>%
  group_by(associated_rodent_record_id) %>%
  slice(which.max(positive)) %>%
  ungroup()

## translate species name into genus/species names

host <- host %>%
  separate(scientificName, into = c("genus", "species", "other"), sep = " ")

## impute negative pathogen



## output file for analysis
combined_data <- list(studies = studies,
                      host = host,
                      pathogen = pathogen)

write_rds(combined_data, file="combined_data_highres.rds")





    