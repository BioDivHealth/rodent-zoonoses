setwd("../data")

## Load packages
if (!require("pacman")) install.packages("pacman")

pkgs <- c(
  "here",
  "tidyverse",
  "lubridate",
  "ggplot2",
  "maps",
  "mapdata"
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

write_rds(combined_data, file="combined_data.rds")
attach(combined_data)

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

## output file for analysis
combined_data <- list(studies = studies,
                      host = host,
                      pathogen = pathogen)

write_rds(combined_data, file="combined_data.rds")



    