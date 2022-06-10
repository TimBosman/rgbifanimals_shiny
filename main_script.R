################################### Main Function ############################################
# Set working directory to directory where the R-script is saved
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # requires installation of package "rstudioapi"
#load in the functions
source("functions.R")
# Create a rastered world
tr <- create_rastered_world("inputs/tr.rdata")
# Read a species list
df <- readRDS("inputs/BOLDigger_Species_Location.rds")
# Read coordinates file
Coordinates <- read.csv("inputs/Coordinates.csv")
# Find for every location the shortest path to the species observation 
sapply(df$Specieslist, function(species){
  print(species)
  tryCatch(find_closest_registered_place(species, Coordinates, tr, "ShortestPath2.csv"), error = function(e)return())
})
 
long <- pivot_longer(df, !Specieslist)
long <- long[long$value > 0, ]

outputfile = "output.csv"

apply(long[1:250,], 1, function(row){
  print(as.character(row[1]))
  if(!file.exists(outputfile)){
    write.clean.csv(c("species", "location","seqsfound", "straightdis", "seadis","pointcalculated"), outputfile)
  }
  samplelocation <- Coordinates[as.character(Coordinates$Observatory.ID) == as.character(row[2]), c("Longitude", "Latitude")]
  tryCatch({occurrence_data <- check_occurrence_data(row[1])},
           error = function(i) {write.clean.csv(c(row, 0, 0, NA), outputfile);return()}
  )
  # Remove duplicates
  occurrence_data <- occurrence_data[!duplicated(occurrence_data),]
  # step1: calculate all distances to samplelocation
  distances <- distm(samplelocation, occurrence_data, fun = distVincentyEllipsoid)
  closest <- occurrence_data[which.min(distances),]
  # step2: Calculate the length through sea for the closest point
  tryCatch({
    path <- shortestPath(tr, sp_format(samplelocation), sp_format(closest), output = "SpatialLines")
    sea_dist <-  geosphere::lengthLine(path)},
    error = function(i) {write.clean.csv(c(row, min(distances), NA, 0), outputfile);return()}
  )
  
  occurrence_data <- occurrence_data[distances <= sea_dist,]
  if(nrow(occurrence_data) <= 1){
    write.clean.csv(c(row,  min(distances), sea_dist, 1), outputfile);return()
  }
  # Save the number of points for which the distance will be calculated
  pointscalculated <- nrow(occurrence_data)
  sea_dists <- sapply(1:nrow(occurrence_data), function(i) {
    path <- shortestPath(tr, sp_format(samplelocation),
                         sp_format(occurrence_data[i,]),
                         output = "SpatialLines")
    return(geosphere::lengthLine(path))
  })
  write.clean.csv(c(row, min(distances), min(sea_dists), pointscalculated), outputfile)
})

 