#------------------------------------------------------------------------------#
# Loading & formatting a rodent phylogeny & calculating phylogenetic distances #
#------------------------------------------------------------------------------#

# 0. Session details ----
# R version 4.3.1 (2023-06-16)
# Running under: Windows 11 x64 (build 22631)

# 1. Load packages ----
pacman::p_load(here,
               stringr,
               dplyr,
               rotl, # for tree of life phylogenies
               ape)  # for phylogenetic analyes

# 2. Load data ----

# Rodent-pathogen data: 
# https://github.com/DidDrog11/scoping_review/tree/main/data_clean
rodent.pathogen.data <- readRDS(here("data", "wide_pathogen.rds")) %>%  # wide format
  as.data.frame() %>%
  sf::st_as_sf()

# 3. Generate a rodent phylogeny ----
