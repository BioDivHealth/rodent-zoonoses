#----------------------------------#
# Explore predictors of prevalence #
#----------------------------------#

# 0. Session details ----
# R version 4.4.0 (2024-04-24)
# Running under: Windows 11 x64 (build 22631)

# Script contents:
# • Load rodent phylogeny, rodent-pathogen & harmonised species names data
# • Filter rodent-pathogen data for pathogen of choice & structure for plotting
# • Calculate phylogenetic distances to reservoir host
# • Plot relative community abundance & phylogenetic distance as prevalence 
#   predictors & save figure

# 1. Load packages ----
pacman::p_load(dplyr, ggplot2, sf, here, stringr, lubridate, ape)

# 2. Load useful custom functions ----
source(here("scripts", "00_useful_functions.R")) #loaded as list `myfuncs`

# 3. Load data ----

# Rodent-pathogen data, from: 
# https://github.com/DidDrog11/scoping_review/tree/main/data_clean
rodent.pathogen.data <- readRDS(here("data", "wide_pathogen.rds")) %>%  # wide format
  as.data.frame() %>%
  sf::st_as_sf()

# Rodent phylogeny, generated in 'scripts/01_rodent_phylogeny.R'
rodent.phylogeny <- read.tree(here("data", "rodentia.mrca.tre"))

# Harmonised rodent species names (matched to names in above phylogeny), generated in 'scripts/01_rodent_phylogeny.R'
harmonised.sp.names <- readRDS(here("data", "harmonised.rodent.sp.names.rds"))

# 4. Format data ----

# Function to filter data by pathogen of choice
select.pathogen <- function(data, pathogen, pathogen.col) {
  # pathogen is a string naming the pathogen of interest
  # pathogen.col is the string name of the pathogen name containing column
  
  # Filter data to have paired abundance-pathogen dataset
  filtered.data <- data %>%
    # Get only trap events where testing for pathogen of choice
    filter(str_detect(!!sym(pathogen.col), pathogen)) %>% 
    # Separate classification into genus and species
    tidyr::separate(col = classification, into = c("genus", "species"), sep = " ", remove = F) %>% 
    # Set species as NA if abbreviated as sp.  
    mutate(species = ifelse(species == "sp.", NA, species)) %>% 
    # Get relevant cols
    select(unique_id, year_trapping, month, country, iso3c, region, town_village,
           classification, genus, species, record_id, name,
           # get all columns containing info on the specified pathogen
           matches(pathogen),
           # get all habitat relevant columns
           contains("habitat")
    )
  
  return(filtered.data)
}


# Function to clean pathogen-specific data & format into hierachical structure
format.hierachy <- function(data) {
  
  data %>%
    # Use strict definition of hierachical grouping:
    filter(
      # Remove rows with missing or non-specific year
      !is.na(year_trapping) & grepl("-", year_trapping) == FALSE &
        # Remove rows with missing month & town/village
        !is.na(month) & !is.na(town_village) & 
        # Remove rows with missing rodent species
        !is.na(species)
    ) %>%
    ### ** Generalise the below **
    # Rename all columns containing 'path_1_tested' to 'abundance'
    rename_with(
      ~ if_else(str_detect(., "path_1_tested"), "abundance", .),
      everything()
    ) %>%
    # Rename all columns containing 'pcr_path_1_positive' to 'pcr.pos'
    rename_with(
      ~ if_else(str_detect(., "pcr_path_1_positive"), "pcr.pos", .),
      everything()
    ) %>% 
    # Rename all columns containing 'culture_path_1_positive' to 'culture.pos'
    rename_with(
      ~ if_else(str_detect(., "culture_path_1_positive"), "culture.pos", .),
      everything()
    ) %>% 
    # Rename all columns containing 'ab_ag_path_1_positive' to 'abag.pos'
    rename_with(
      ~ if_else(str_detect(., "ab_ag_path_1_positive"), "abag.pos", .),
      everything()
    ) %>% 
    # Rename all columns containing 'histo_path_1_positive' to 'histo.pos'
    rename_with(
      ~ if_else(str_detect(., "histo_path_1_positive"), "histo.pos", .),
      everything()
    ) %>% 
    # Calculate relative community abundance for each trapping event
    group_by(unique_id, #study
             year_trapping, month, #sampling timepoint for given study
             town_village #site
    ) %>%
    mutate(
      # Total abundance at each site
      site.tot.abund = sum(abundance),
      # Site-level ID
      site.id = paste(unique_id, year_trapping, month, town_village, row_number(), sep = "_") ) %>%
    # Relative abundance of each species at each site (group by host species)
    group_by(site.id, classification) %>%
    mutate(
      site.species.abund = sum(abundance),
      site.species.rel.abun = site.species.abund / site.tot.abund,
      # Calculate prevalence for different methods -- ** generalise this **
      prevalence.pcr   = pcr.pos / abundance,
      prevalence.cult  = culture.pos / abundance,
      prevalence.abag  = abag.pos / abundance#,
      # prevalence.histo = histo.pos / abundance
    ) %>% 
    # Remove old columns
    select(!contains(c("tested", "positive")))
}

# Run pathogen filtering function
path.data <- select.pathogen(data = rodent.pathogen.data, pathogen = "lassa", pathogen.col = "name")

# Run pathogen data cleaning & hierarchy formatting function
path.data.clean <- format.hierachy(data = path.data)

# Add harmonised rodent species names to data (matching nomenclature of phylogeny)
path.data.clean <- path.data.clean %>% 
  left_join(., harmonised.sp.names %>% select(classification, tree.names)) %>% 
  # Remove any rows with no match to species in phylogeny
  filter(!is.na(tree.names)) %>%
  # Remove species not present in the rodent tree
  filter(tree.names %in% rodent.phylogeny$tip.label)

# Are all rodent species present in the phylogeny?
# all(path.data.clean$tree.names %in% rodent.phylogeny$tip.label)

# Get a pairwise (relative) distance matrix
dm <- myfuncs$get.phydist.mat(phylogeny = rodent.phylogeny, rel.dist = TRUE)

# Estimate phylogenetic distances from each species to the reservoir host
reservoir.host <- "Mastomys_natalensis"

path.data.clean <- path.data.clean %>% 
  mutate(phy.dist.reservoir = myfuncs$get.phydist(distance.matrix = dm, 
                                                   sp1 = reservoir.host, 
                                                   sp2 = tree.names))

# 5. Plot community abundance vs. lassa prevalence ----
path.data.clean %>% 
  # Put data in long format & make use of all prevalence data (via all assays)
  tidyr::pivot_longer(cols = c(prevalence.pcr, prevalence.cult, prevalence.abag), 
                      names_to = "assay", values_to = "prevalence") %>% 
  mutate(
    # Rename assay to remove 'prevalence' prefix
    assay = str_replace(assay, "prevalence\\.", ""),
    # Tidy up assay naming
    assay = case_when(
      assay == "pcr"  ~ toupper(assay),
      assay == "abag" ~ "Ab/Ag",
      assay == "cult" ~ "Culture"
    )) %>% 
  # Make plot
  ggplot() +
  geom_point(aes(site.species.rel.abun, prevalence * 100, 
                 # Colour points by phylogenetic distance to reservoir
                 col = phy.dist.reservoir, 
                 # Size points by total site abundance
                 size = site.tot.abund,
                 # Shape points by assay used to detect pathogen
                 shape = assay), 
             alpha = 0.5) + 
  theme_minimal() + 
  labs(title = "Relative community abundance & phylogenetic distance as predictors of lassa prevalence",
       x     = "Relative community abundance of rodent species",
       y     = "Prevalence of lassa infection (%)",
       col   = "Phylogenetic distance to reservoir",
       shape = "Assay",
       size  = "Total community abundance") +
  # Use viridis colour scale
  scale_color_viridis_c() + 
  # Format text size & face
  theme(
    text          = element_text(size = 10, family = "serif"), 
    plot.title    = element_text(size = 10, family = "serif", face = "bold"),
    plot.subtitle = element_text(size = 10, family = "serif", face = "italic")
  )

# Save plot
ggsave("figures/lassa_relabundance.pdf", width = 6, height = 5, units = "in")
