## Load packages
if (!require("pacman")) install.packages("pacman")

pkgs <- c(
  "here",
  "tidyverse",
  "lubridate",
  "dplyr"
)

source(here("scripts", "00_useful_functions.R")) #loaded as list `myfuncs`

pacman::p_load(pkgs, character.only = T)

## if multiple in put sheets, load and merge sheets
##harry
studies_h = read.csv("./data/aren_hant_data/studies_h.csv")
studies_h = studies_h[0:8]
host_h = read.csv("./data/aren_hant_data/host_h.csv")
host_h = host_h[0:13]
pathogen_h = read.csv("./data/aren_hant_data/pathogen_h.csv")
pathogen_h = pathogen_h[0:12]
##david
studies_d = read.csv("./data/aren_hant_data/studies_d.csv")
host_d = read.csv("./data/aren_hant_data/host_d.csv")
pathogen_d = read.csv("./data/aren_hant_data/pathogen_d.csv")
pathogen_d = pathogen_d[0:12]

##ana
studies_a = read.csv("./data/aren_hant_data/studies_a.csv")
studies_a = studies_a[0:8]
host_a = read.csv("./data/aren_hant_data/host_a.csv")
host_a = host_a[0:13]
pathogen_a = read.csv("./data/aren_hant_data/pathogen_a.csv")


studies <-  rbind(studies_h, studies_d, studies_a)
host <-  rbind(host_h, host_d, host_a)
pathogen <-  rbind(pathogen_h, pathogen_d, pathogen_a)


## Filter for site resolution

studies <- subset(studies, data_resolution == "site-session")

# keep only host and pathogen data with matching study ID from this subset

host <- host[host$study_id %in% studies$study_id, ]
pathogen <- pathogen[pathogen$study_id %in% studies$study_id, ]

# Keep only 1 pathogen row per rodent
## check for multiple pathogen entries for each rodent_id and take max positives

pathogen <- pathogen %>%
  group_by(associated_rodent_record_id) %>%
  slice(which.max(positive)) %>%
  ungroup()

## Impute pathogen negatives
## get rodent record IDs with no corresponding pathogen entry
host_impute <- host[!(host$rodent_record_id %in% pathogen$associated_rodent_record_id), ]

## Get relevant columns from this dataframe
new_pathogen <- data.frame(host_impute$rodent_record_id,
                           host_impute$study_id,
                           host_impute$scientificName,
                           host_impute$individualCount)

# rename columns to match pathogen sheet
new_pathogen <- new_pathogen %>% 
  rename(
    associated_rodent_record_id = host_impute.rodent_record_id,
    study_id = host_impute.study_id,
    associatedTaxa = host_impute.scientificName,
    tested = host_impute.individualCount
  )

new_pathogen <- merge(new_pathogen, pathogen[, c("study_id", "family", "scientificName", "assay")], by = "study_id", all.x = TRUE)

new_pathogen <- unique(new_pathogen) ### THINK ABOUT NEGATIVE ARENAVIRUS/HANTAVIRUS TESTS WHERE BOTH ARE TESTED BUT ONLY 1 REPORTED


## group individual-level data to (roughly) monthly data

## subset individual level data


host_individual <- subset(host, study_id == "h_129")
host_individual$eventDate <- as.Date(host_individual$eventDate)
host_individual <- host_individual %>% 
  group_by(locality) %>%
  group_modify(~ group_time_series(.x)) %>%
  mutate(group = paste(locality, group, sep = "_")) %>%
  ungroup()

## merge with pathogen data

pathogen <- pathogen %>% 
  rename(
    rodent_record_id = associated_rodent_record_id
  )

host_individual <- merge(host_individual,pathogen,by="rodent_record_id")

## group and summarise to obtain summarised data

host_individual <- host_individual %>%
  group_by(group, scientificName.x) %>%
  reframe(rodent_record_id = first(rodent_record_id),
            study_id.x = first(study_id.x),
            start_date = min(eventDate),
            end_date = max(eventDate),
            locality = first(locality),
            country = first(country),
            verbatimLocality = first(verbatimLocality),
            coordinate_resolution = first(coordinate_resolution),
            decimalLatitude = first(decimalLatitude),
            decimalLongitude = first(decimalLongitude),
            individualCount = sum(as.numeric(individualCount)),
            trapEffort = first(trapEffort),
            trapEffortResolution = first(trapEffortResolution),
            pathogen_record_id = first(pathogen_record_id),
            study_id.y = first(study_id.y),
            associatedTaxa = first(associatedTaxa),
            family = first(family),
            scientificName.y = first(scientificName.y),
            assay = first(assay), #### this is not accurate with current code ###
            tested = sum(as.numeric(tested)),
            positive = sum(as.numeric(positive)),
            negative = sum(as.numeric(negative)),
            number_inconclusive = first(number_inconclusive),
            note = first(note)
            )

host_individual <- host_individual[, -1]

## translate summarised host dataset date into start/end columns

host_summarised <- subset(host, study_id != "h_129")

host_summarised <- host_summarised %>%
  separate(eventDate, into = c("start_date", "end_date"), sep = "/")

host_monthly <- host_summarised %>%
  filter(str_detect(start_date, "^\\d{4}-\\d{2}") | str_detect(start_date, "^\\d{4}-\\d{2}-\\d{2}"))

# Add "-01" to strings in "yyyy-mm" format
host_monthly$start_date <- ifelse(grepl("^\\d{4}-\\d{2}$", host_monthly$start_date), 
                                  paste0(host_monthly$start_date, "-01"), 
                                  host_monthly$start_date)

host_monthly$end_date <- ifelse(grepl("^\\d{4}-\\d{2}$", host_monthly$end_date), 
                                paste0(host_monthly$end_date, "-01"), 
                                host_monthly$end_date)

## Merge pathogen and host summarised data

host_monthly <- merge(host_monthly,pathogen,by="rodent_record_id")

## rbind individual and summarised data

host_monthly <- rbind(host_monthly, host_individual)

## convert all dates to yyyy-mm-dd and convert to date object

host_monthly$start_date <- as.Date(host_monthly$start_date)
host_monthly$end_date <- as.Date(host_monthly$end_date)

## Calculate duration of trapping
host_monthly$period <- host_monthly$end_date - host_monthly$start_date
host_monthly$period[is.na(host_monthly$period)] <- 0

## filter for less than 2 month trapping time
host <- host_monthly %>% filter(period <= 60)

## filter for multiple species
host <- host %>%
  group_by(study_id.x) %>%
  filter(n_distinct(scientificName.x) > 1) %>%
  ungroup()

## convert Ana's coordinates to numeric

host$decimalLatitude <- as.numeric(host$decimalLatitude)
host$decimalLongitude <- as.numeric(host$decimalLongitude)



## impute negative pathogen sheet entries



## if( %in% studies_to_impute$study_id) {



## translate species name into genus/species names

host_path_wide <- host %>%
  separate(scientificName.x, into = c("genus", "species", "other"), sep = " ")

# save to rds for phylogeny

write_rds(host_path_wide, file="./data/host_path_wide.rds")


### how many communities do we have?

unique(host[c('start_date','locality')])
