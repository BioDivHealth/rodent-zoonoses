setwd("../data")
if (!require("pacman")) install.packages("pacman")

pkgs <- c(
  "here",
  "tidyverse",
  "lubridate"
)


pacman::p_load(pkgs, character.only = T)

studies = read.csv("studies.csv")
host = read.csv("host.csv")
pathogen = read.csv("pathogen.csv")


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



