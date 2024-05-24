#-----------------------------------------#
# Loading & formatting a rodent phylogeny #
#-----------------------------------------#

# 0. Session details ----
# R version 4.4.0 (2024-04-24)
# Running under: Windows 11 x64 (build 22631)

# Script contents:
# • Load mammal phylogeny, prune to get most recent common ancestor & save as 
#   .tre file for subsequent analyses
# • Harmonise rodent species names between rodent-pathogen data and tree & save
#   as .rds file
# • Plot mammal phylogeny & show position of rodent species in the data, save
#   as figure

# 1. Load packages ----
pacman::p_load(here,
               stringr,
               dplyr,
               ape,    # for phylogenetic analyes
               ggtree, # for plotting trees
               ggplot2)  

# 2. Load data ----

# Rodent-pathogen data (for testing phylogenetic workflow)
# https://github.com/DidDrog11/scoping_review/tree/main/data_clean
rodent.pathogen.data <- readRDS(here("data", "wide_pathogen.rds")) %>%  # wide format
  as.data.frame() %>% 
  sf::st_as_sf() #to avoid error: Column `geometry` is a `sfc_POINT/sfc` object

# Load mammal phylogeny: https://data.vertlife.org -> mammaltree -> Completed_5911sp_topoCons_FBDasZhouEtAl.zip
# I'm using one tree of 1,000 bootstrap replicates; can extend this once workflow established
mammal.tree <- read.tree(here("data", "MamPhy_BDvr_Completed_5911sp_topoCons_FBDasZhouEtAl_v2_tree0000.tre"))

# Summarise the tree
str(mammal.tree)

# 3. Generate rodent phylogeny ----

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
rodent.species.corrected <- rodent.species %>% 
  mutate(
    tree.names = case_when(
      # Acomys_chudeaui = Acomys chudeaui, proposed to be a junior synonym for Acomys airensis: https://academic.oup.com/biolinnean/article/98/1/29/2235991
      tree.names == "Acomys_chudeaui" ~ "Acomys_airensis",
      # theresae was misspelt
      tree.names == "Crocidura_thereseae" ~ "Crocidura_theresae", 
      # Dipodillus is sometimes classified as a subgenus of Gerbillus
      tree.names == "Dipodillus_campestris" ~ "Gerbillus_campestris",
      # guineae was misspelt 
      tree.names == "Gerbilliscus_guinea" ~ "Gerbilliscus_guineae",
      # kempi was misspelt
      tree.names == "Gerbilliscus_kempii" ~ "Gerbilliscus_kempi",
      # Hylomyscus_simus: not found anywhere online, so set name to NA
      tree.names == "Hylomyscus_simus" ~ NA,
      # Keep all other names the same
      TRUE ~ tree.names
    ),
  ) %>% 
  # Drop species with NA tree names
  filter(!is.na(tree.names)) %>% 
  # Keep only unique species names according to their names in the tree
  filter(duplicated(tree.names) == FALSE)

# Ensure all species are now in the tree
rodent.species.corrected$in.tree <- ifelse(rodent.species.corrected$tree.names %in% mammal.tree$tip.label,
                                           TRUE, FALSE)

sum(rodent.species.corrected$in.tree == FALSE)  # now 0 missing from tree

# Save rodent species tree names for use in other scripts
saveRDS(rodent.species.corrected, here("data", "harmonised.rodent.sp.names.rds"))

# Find the node number for the clade containing all species in the dataset
rodentia.node <- getMRCA(mammal.tree, tip = rodent.species.corrected$tree.names)

# Extract the clade using the node number 
rodentia.clade <- extract.clade(mammal.tree, rodentia.node)

# Summarise the extracted clade
print(rodentia.clade)

# Save the rodent phylogeny
write.tree(rodentia.clade, here("data", "rodentia.mrca.tre"))

# 4. Plot phylogeny ----

# Create variable to colour tips if species present in the data
rodentia.clade$in.data <- ifelse(rodentia.clade$tip.label %in% rodent.species.corrected$tree.names, 
       "purple", NA)
rodentia.clade$in.data <- as.factor(rodentia.clade$in.data)

# Plot the tree using the ggtree package
rodentia.clade %>% 
  ggtree() + 
  # Add an evolutionary distance scale bar
  geom_treescale(x = 0, y = 45, col = "midnightblue") + 
  # Add points for species in the dataset (NB: doesnt work without data$column syntax)
  geom_tippoint(color = rodentia.clade$in.data)
