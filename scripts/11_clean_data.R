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

## convert all dates to yyyy-mm-dd and convert to date object









