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

# 1.5. Load useful custom functions ----
source(here("scripts", "00_useful_functions.R")) #loaded as list `myfuncs`

# 2. Load data ----

# Rodent-pathogen data (for testing phylogenetic workflow)
# https://github.com/DidDrog11/scoping_review/tree/main/data_clean
host.pathogen.data <- readRDS(here("data", "host_path_wide.rds")) %>%  # wide format
  as.data.frame()

# Load mammal phylogeny: https://data.vertlife.org -> mammaltree -> Completed_5911sp_topoCons_FBDasZhouEtAl.zip
# I'm using one tree of 1,000 bootstrap replicates; can extend this once workflow established
mammal.tree <- read.tree(here("data", "./phylogenetics/MamPhy_BDvr_Completed_5911sp_topoCons_FBDasZhouEtAl_v2_tree0000.tre"))

# Load species taxonomical rankings (https://github.com/n8upham/MamPhy_v1/blob/master/_DATA/taxonomy_mamPhy_5911species_toPublish.csv)
taxa <- read.csv(here("data", "./phylogenetics/taxonomy_mamPhy_5911species_toPublish.csv")) 

# Summarise the tree
str(mammal.tree)

# 3. Generate rodent phylogeny ----

# Format data to get host names in same format as trees
host.species <- host.pathogen.data %>%
  mutate(
    # Set species as NA if abbreviated as sp.  
    species = ifelse(species == "sp.", NA, species),
    species = ifelse(species == "spp.", NA, species),
    # Capitalise first letter of genus
    genus = paste0(toupper(substr(genus, 1, 1)), substr(genus, 2, nchar(genus)))
  ) %>% 
  # Remove NA species
  filter(!is.na(species)) %>% 
  # Create classification that matches tree tip labels in format Genus_species
  mutate(tree.names = paste(genus, species, sep = "_")) %>% 
  # Return select columns
  select(associatedTaxa,tree.names, genus, species) %>% 
  # Keep only unique species names
  distinct() %>% 
  # Order alphabetically
  arrange(associatedTaxa)

# Which of our rodent species are in the tree?
host.species$in.tree <- ifelse(host.species$tree.names %in% mammal.tree$tip.label,
                                 TRUE, FALSE)

sum(host.species$in.tree == FALSE)  #3 missing from tree

# For those missing from the tree, which same-genus species does the tree contain?
missing.species <- host.species %>% 
  filter(in.tree == FALSE)

for(i in 1:length(missing.species)){
  same.genus.species <- mammal.tree$tip.label[str_detect(mammal.tree$tip.label, missing.species$genus[i])] %>% 
    sort() # sort output alphabetically
  print(paste0("Missing species in rodent data: ", missing.species$tree.names[i]))
  print("Tree species in the same genus:")
  print(same.genus.species)
}

# Correct names in the rodent data after manual comparison with tree names
host.species.corrected <- host.species %>% 
  mutate(
    tree.names = case_when(
      # azarae was misspelt
      tree.names == "Akodon_azare" ~ "Akodon_azarae",
      # Clethrionomys was misspelt and blank vole is also referred to as myodes
      tree.names == "Clethryonomys_glareolus" ~ "Myodes_glareolus", 
      # Icttidomys mexicanus is sometimes classed as Spermophilus
      tree.names == "Spermophilus_mexicanus" ~ "Ictidomys_mexicanus",
      # Keep all other names the same
      TRUE ~ tree.names
    ),
  ) %>% 
  # Drop species with NA tree names
  filter(!is.na(tree.names)) %>% 
  # Keep only unique species names according to their names in the tree
  filter(duplicated(tree.names) == FALSE)

# Ensure all species are now in the tree
host.species.corrected$in.tree <- ifelse(host.species.corrected$tree.names %in% mammal.tree$tip.label,
                                           TRUE, FALSE)

sum(host.species.corrected$in.tree == FALSE)  # now 0 missing from tree

# Save rodent species tree names for use in other scripts
if (!(file.exists(here("data", "./phylogenetics/harmonised.host.sp.names.rds")))) {
  saveRDS(host.species.corrected, here("data", "./phylogenetics/harmonised.host.sp.names.rds"))
}

# Find the node number for the clade containing all species in the dataset
node.in.data <- getMRCA(mammal.tree, tip = host.species.corrected$tree.names)

# Extract the clade using the node number
pruned.mammal.tree <- extract.clade(mammal.tree, node.in.data)

# Remove species in the data that are not within our clade of interest

# Get species name from tree
mammal.spp <- data.frame(Species_Name = pruned.mammal.tree$tip.label) %>% 
  # Join to taxonomic family & order data
  left_join(taxa %>% select(Species_Name, fam, ord), by = "Species_Name")
  # Get only species in the order Rodentia

# Filter host species in the data to include only those in our dataset's clade
pruned.host.spp <- host.species.corrected %>% 
  filter(tree.names %in% mammal.spp$Species_Name)

# nrow(rodent.species.corrected)  #before pruning (include non-rodents)
# nrow(pruned.host.spp)           #after pruning (only rodents)

# Find the node number for the clade containing all our host species in the dataset
host.node <- getMRCA(mammal.tree, tip = pruned.host.spp$tree.names)

# Extract the clade using the node number 
host.clade <- extract.clade(mammal.tree, host.node)

# Summarise the extracted clade
str(host.clade)

# Save the rodent phylogeny
if (!(file.exists(here("data", "./phylogenetics/host.mrca.tree")))) {
  write.tree(host.clade, here("data", "./phylogenetics/host.mrca.tree"))
}

# Remove large objects from the environment
rm(mammal.tree)

# 4. Plot phylogeny ----

# Create variable to colour tips if species present in the data
host.clade$in.data <- ifelse(host.clade$tip.label %in% pruned.host.spp$tree.names, 
       "purple", NA)
host.clade$in.data <- as.factor(host.clade$in.data)

# Plot the phylogeny
ggtree(host.clade, size = 0.05) + 
  geom_tippoint(color = host.clade$in.data, size = 3) + 
  geom_treescale(x = 0, y = 2200, col = "midnightblue")

# Save the phylogeny
ggsave(here("figures", "host.phylogeny.pdf"), width = 6, height = 6, units = "in")

### get distance matrix

# Get a pairwise (relative) distance matrix
dm <- myfuncs$get.phydist.mat(phylogeny = pruned.mammal.tree, rel.dist = TRUE)

species <- host.species.corrected$tree.names

## check all the names are present in the matrix
missing_names <- setdiff(species, rownames(dm))
if(length(missing_names) > 0) {
  warning("The following names are not in the distance matrix and will be ignored: ", paste(missing_names, collapse = ", "))
}

# only keep our species of interest
filtered.dm <- dm[species, species]

# save final matrix
write_rds(filtered.dm, file="./phylogenetics/phylo_dist_matrix.rds")
