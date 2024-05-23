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
               ape)  # for phylogenetic analyes

# 2. Load data ----

# Rodent-pathogen data (for testing phylogenetic workflow)
# https://github.com/DidDrog11/scoping_review/tree/main/data_clean
rodent.pathogen.data <- readRDS(here("data", "wide_pathogen.rds")) %>%  # wide format
  as.data.frame() %>% 
  sf::st_as_sf() #to avoid error: Column `geometry` is a `sfc_POINT/sfc` object

# 3. Generate a rodent phylogeny ----

# https://data.vertlife.org -> mammaltree -> Completed_5911sp_topoCons_FBDasZhouEtAl.zip

# I'm using one tree of 1,000 bootstrap replicates; can extend this once workflow established
mammal.tree <- read.tree(here("data", "MamPhy_BDvr_Completed_5911sp_topoCons_FBDasZhouEtAl_v2_tree0000.tre"))

# Summarise the tree
str(mammal.tree)
# plot.phylo(mammal.tree, cex = 0.2)  #NB: runs slowly

# Format rodent data to get host names in same format as tree's
rodent.species <- rodent.pathogen.data %>%
  # Temporarily drop the 'geometry' column
  sf::st_drop_geometry() %>%
  # Separate host classification into genus and species
  tidyr::separate(col = classification, into = c("genus", "species"), sep = " ", remove = F) %>% 
  mutate(
    # Set species as NA if abbreviated as sp.  
    species = ifelse(species == "sp.", NA, species),
    # Capitalise first letter of genus
    genus = paste0(toupper(substr(genus, 1, 1)), substr(genus, 2, nchar(genus)))
  ) %>% 
  # Remove NA species
  filter(!is.na(species)) %>% 
  # Create classification that matches tree tip labels
  mutate(tree.names = paste(genus, species, sep = "_")) %>% 
  # Return select columns
  select(classification, tree.names, genus, species) %>% 
  # Only unique species names
  distinct() %>% 
  # Order alphabetically
  arrange(classification)

# Which of our rodent species are in the tree?
rodent.species$in.tree <- ifelse(rodent.species$tree.names %in% mammal.tree$tip.label,
                                 TRUE, FALSE)

sum(rodent.species$in.tree == FALSE)  #6 missing from tree

# For those missing from the tree, which same-genus species does the tree contain?
missing.species <- rodent.species %>% 
  filter(in.tree == FALSE)

for(i in 1:length(missing.species)){
  same.genus.species <- mammal.tree$tip.label[str_detect(mammal.tree$tip.label, missing.species$genus[i])] %>% 
    sort() # sort output alphabetically
  print(paste0("Missing species in rodent data: ", missing.species$tree.names[i]))
  print("Tree species in the same genus:")
  print(same.genus.species)
}

# Correct names in the rodent data after manual comparison with tree names
# ** TO DO: search for the 2 remaining missing species **
rodent.species.corrected <- rodent.species %>% 
  mutate(
    tree.names = case_when(
      # theresae was misspelt
      tree.names == "Crocidura_thereseae" ~ "Crocidura_theresae", 
      # Dipodillus is sometimes classified as a subgenus of Gerbillus
      tree.names == "Dipodillus_campestris" ~ "Gerbillus_campestris",
      # guineae was misspelt 
      tree.names == "Gerbilliscus_guinea" ~ "Gerbilliscus_guineae",
      # kempi was misspelt
      tree.names == "Gerbilliscus_kempii" ~ "Gerbilliscus_kempi", 
      # Keep all other names the same
      TRUE ~ tree.names
    )
  )

# How many non-matching species after correcting names?
num.missing <- ifelse(rodent.species.corrected$tree.names %in% mammal.tree$tip.label,
                      TRUE, FALSE)
sum(num.missing == FALSE)

# Find the node number for the clade containing all species in the dataset
# NB: node of MRCA for all rodents = 2392; from https://vertlife.org/data/mammals/
rodentia.node <- getMRCA(mammal.tree, 
                         # ** to do: ensure all species are matched to the tree then change the below **
                         tip = rodent.species.corrected$tree.names[rodent.species.corrected$in.tree])

# Extract the clade using the node number 
rodentia.clade <- extract.clade(mammal.tree, rodentia.node)

# Summarise the extracted clade
print(rodentia.clade)
plot(rodentia.clade, cex = 0.2)  # NB: runs slowly
