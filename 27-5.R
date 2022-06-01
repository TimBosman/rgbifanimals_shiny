# Load packages
#library("raster")
library("gdistance")
library("maptools")
library("dplyr")
# library("rgdal")
# library("jsonlite")
require("geosphere")
# require("curl")
require("rgbif")
# library("maps")
library("ggplot2")
library("ggrepel")
library("tidyr")
# library("units")
# library("sf")
# library("ggspatial")
#require(parallel)

create_rastered_world <- function(filename){
  # if the file is already made, it will load the file and return it
  if(file.exists(filename)){
    load(filename)
    return(tr)
  }
  # Load a map
  data(wrld_simpl) #use wrld_simpl from the maptools package
  
  # Generate a scaffold for the raster file
  world_crs <- crs(wrld_simpl)
  worldshp <- spTransform(wrld_simpl, world_crs)
  ras <- raster(nrow=1200, ncol=1200)
  
  # Generate a raster file
  worldmask <- rasterize(worldshp, ras)
  worldras <- is.na(worldmask) # inverse water and land, so ocean becomes 1 and land 0
  worldras[worldras==0] <- 999 # set land to 999
  
  # Create a Transition object from the raster
  tr <- transition(worldras, function(x) 1/mean(x), 16)
  tr = geoCorrection(tr, scl=FALSE)
  save(tr, file = filename)
  return(tr)
}

get_occurence_data <- function(species){
  res = occ_data(scientificName = species, hasCoordinate = TRUE, limit = 200000)
  res <- res$data[, c('decimalLongitude', 'decimalLatitude')]
  #rename the column names
  colnames(res) <- c('Longitude', 'Latitude')
  # Remove occurances where longitude or latitude is NA
  res <- res[!is.na(res$Latitude) & !is.na(res$Longitude),]
  return(res)
}

check_occurence_data <- function(species){
  filename = paste0("OccuranceData/", species, ".csv")
  if(file.exists(filename)){
    res <- read.csv(filename, header = TRUE)
    return(res)
  } else {
    res <- get_occurence_data(species)
    write.csv(res, file = filename, quote = FALSE, row.names = FALSE)
    return(res)
  }
}

plot_distribution <- function(Coordinates, res, plotname, title){
  ggplot() +
    geom_polygon(aes(x = long, y = lat, group = group), data = map_data("world")) +
    geom_hex(aes(x= Longitude, y = Latitude), data = res) +
    # coord_cartesian(xlim = c(-30, 60), ylim = c(30,70)) +
    # geom_label_repel(aes(x= Longitude, y = Latitude, label = Observatory.ID), data = Coordinates) +
    geom_point(aes(x= Longitude, y = Latitude), data = Coordinates, col = "red") + 
    ggtitle(title)
  ggsave(plotname)
}

find_closest_registered_place <- function(species, Coordinates, tr, outputfile, plot=TRUE){
  # Check the occurence data, If there is an error there it catches it and writes an error in the output file
  tryCatch(res <- check_occurence_data(species),
           error = function(errormessage){
             write.table(paste(c(species, rep("", 15), "ERROR while getting the occurence data"), 
                               collapse = ","), file = outputfile, append = TRUE, quote = FALSE, 
                         col.names = FALSE, row.names = FALSE)
           })
  obs = nrow(res)
  # Plot the distribution
  plotname = paste0("plots/", species, ".jpeg")
  if(plot & !file.exists(plotname)){
    plot_distribution(Coordinates, res, plotname, species)
  }
  # Remove duplicate coordinates
  res <- res[!duplicated(res),]
  uobs = nrow(res)
  # Calculate the distances between all the points and all the locations
  distances <- distm(res[, c("Longitude", "Latitude")], Coordinates[, c("Longitude", "Latitude")], 
                     fun = distVincentyEllipsoid)
  # Find the closest point to every sampling location
  shortest <- round(apply(distances, 2, function(dist) min(dist)), 0)
  # Add headers to the output file if it is not made already
  if(!file.exists(outputfile)){
    write.table(paste(c("Speciesname", "Total observations", "Unique locations", 
                as.character(Coordinates$Observatory.ID), "Errors"), collapse = ","), 
                file = outputfile, quote = FALSE, col.names = FALSE, row.names = FALSE)
  }
  # Writing results to the file
  write.table(paste(c(species, obs, uobs, shortest, ""), collapse = ","), file = outputfile, append = TRUE, 
              quote = FALSE, col.names = FALSE, row.names = FALSE)
}

find_shortest_route_in_sea <- function(location, observations, tr){
  #Remove duplicates
  observations <- observations[!duplicated(observations),]
  # Add column with distance
  if(nrow(observations) > 10){
    observations$distance <- pmax(abs(observations$Longitude - location$Longitude), 
                         abs(observations$Latitude - location$Latitude))
    # Get the 5th smallest distance and round it up
    dist = ceiling(sort(observations$distance)[10])
    observations <- observations[observations$distance < dist,]
  }
  # I'm not completely sure about this part of the code. I think the bug should be fixed not avoided :(
  #Change any entries which are considered to be in the same "box" on the map as the sampling location. (observationsolution is 0.3, 0.15), as this would yield an error in the shortestPath function. 
  observations$Latitude <- ifelse(between(observations$Longitude, location$Longitude - 0.3, location$Longitude + 0.3) &
                                  between(observations$Latitude, location$Latitude - 0.15, location$Latitude), 
                                observations$Latitude - 0.16,
                                observations$Latitude)
  observations$Latitude <- ifelse(between(observations$Longitude, location$Longitude - 0.3, location$Longitude + 0.3) &
                                  between(observations$Latitude, location$Latitude, location$Latitude + 0.15), 
                                observations$Latitude + 0.16,
                                observations$Latitude)
  # I don't see why we do this, but okay
  x <- rep("NA", nrow(observations))
  # find the shortest route to every point through the sea
  x <- sapply(1:nrow(observations), function(i) {
    gbif_occurrence <- structure(c(observations$Longitude[i], observations$Latitude[i]), .Dim = 1:2)
    path <- shortestPath(tr, structure(as.numeric(location), .Dim = 1:2), gbif_occurrence, output = "SpatialLines")
    x[i] <- geosphere::lengthLine(path)
  })
  # Find the closest location the point of sampling
  distance <- min(as.numeric(x), na.rm=T)/1000
  return(distance)
}

################################### Main Function ###########################################

# Set working directory to directory where the R-script is saved
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # requires installation of package "rstudioapi"
#Create a rastered world
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
 
# df <- pivot_longer(df, !Specieslist)
# df <- df[df$value > 0, ]
# 
# apply(df[1:2,], 1, function(row){
#   samplelocation <- Coordinates[Coordinates$Observatory.ID == row[2], c("Longitude", "Latitude")]
#   occurence_data <- check_occurence_data(row[1])
#   dist <- find_shortest_route_in_sea(samplelocation, occurence_data, tr)
#   # return(row)
#   return(dist)
# })

