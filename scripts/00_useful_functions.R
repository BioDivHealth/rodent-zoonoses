#------------------#
# Useful functions #
#------------------##

# 0. Session details ----
# R version 4.4.0 (2024-04-24)
# Running under: Windows 11 x64 (build 22631)

# 1. Define a list of useful functions to be sourced in other scripts ----

myfuncs <- list(
  
  # Derive phylogenetic distance matrix given some phylogeny
  get.phydist.mat = function(phylogeny, rel.dist){
    
    # •'phylogeny': a class phylo object with tip labels & edge (branch) lengths
    # • 'rel.dist': a logical, giving a matrix of absolute distances (FALSE) or 
    #               relative distances (TRUE) on the interval [0,1]
    
    # Derive a pairwise phylogenetic distance matrix using ape::cophenetic.phylo()
    distance.matrix <- ape::cophenetic.phylo(phylogeny)
    
    # Transform distance matrix to a relative distance matrix with max value of 1
    if (rel.dist) distance.matrix <- distance.matrix / max(distance.matrix)
    
    return(distance.matrix)
    
  },
  
  
  # Estimate phylogenetic distance between two species
  get.phydist = function(distance.matrix, sp1, sp2){
    
    # • 'distance.matrix': numeric distance matrix generated using ape::cophenetic.phylo()
    # • 'sp1': character string in the form Genus_species
    # • 'sp2': character string in the form Genus_species
    
    # Error handling: check if species name present in distance matrix
    species.names <- colnames(distance.matrix)
    
    if (all(sp1 %in% species.names == FALSE)) {
      stop("Species 1 not named in the input distance matrix")
    }
    
    if (all(sp2 %in% species.names == FALSE)) {
      stop("Species 2 not named in the input distance matrix")
    }  
    
    # Extract & return pairwise phylogenetic distance
    phy.dist <- distance.matrix[sp1, sp2]
    return(phy.dist)
    
  },
  


  glc.files <- function(data, column) {
    # This function returns the glc files you need for your coordinate data
    # Check if the specified column exists in the dataset
    if (!column %in% colnames(data)) {
      stop("Specified column does not exist in the dataset")
    }
    
    # Extract the values from the specified column
    values <- data[[column]]
    
    # Create a sequence of bin floors from -180 to 175 (since 180 would start a new bin)
    bins <- seq(-180, 175, by = 5)
    
    # Find the bin floor for each value
    bin_floors <- sapply(values, function(x) {
      bin <- max(bins[bins <= x])
      return(as.character(bin))
    })
    
    # Remove duplicate bins, sort them, and return as a list of strings
    unique_bin_floors <- unique(bin_floors)
    sorted_bin_floors <- sort(as.numeric(unique_bin_floors))
    
    # Convert sorted bin floors back to character strings
    sorted_bin_floors_strings <- as.character(sorted_bin_floors)
    
    return(sorted_bin_floors_strings)
  }

)