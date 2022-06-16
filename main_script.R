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
  tryCatch(find_closest_registered_place(species, Coordinates, tr, "ShortestPath.csv"), error = function(e)return())
})
 
long <- pivot_longer(df, !Specieslist)
long <- long[long$value > 0, ]

apply(long, 1, function(row){
  samplelocation <- Coordinates[Coordinates$Observatory.ID == row[2], c("Longitude", "Latitude")]
  print(row[1])
  tryCatch({occurrence_data <- check_occurrence_data(row[1])
            find_shortest_route_in_sea(samplelocation, occurrence_data, tr, row, "accurate.csv")},
           error = function(errormessage){
             write.table(paste(c(row, 0, 0, NA), collapse = ","), file = "accurate.csv", 
                         append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE)
           })
})


# results <- read.csv("accurate.csv")
# results <- results[!is.na(results$distance),]
# write.csv(results, "accurate.csv", quote = F, row.names = F)
