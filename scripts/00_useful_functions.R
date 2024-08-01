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
  
  
  get_dec_box <- function(lon, lat) {
    upper_left_lat <- ceiling(lat / 5) * 5
    upper_left_lon <- floor(lon / 5) * 5
    return(c(upper_left_lat, upper_left_lon))
  }

)