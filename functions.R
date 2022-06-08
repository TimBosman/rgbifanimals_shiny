################################### Load Packages ############################################
library("gdistance")
library("maptools")
library("dplyr")
require("geosphere")
require("rgbif")
library("ggplot2")
library("ggrepel")
library("tidyr")

################################### Functions ################################################

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

get_occurrence_data <- function(species){
  res = occ_data(scientificName = species, hasCoordinate = TRUE, limit = 200000)
  res <- res$data[, c('decimalLongitude', 'decimalLatitude')]
  #rename the column names
  colnames(res) <- c('Longitude', 'Latitude')
  # Remove occurrences where longitude or latitude is NA
  res <- res[!is.na(res$Latitude) & !is.na(res$Longitude),]
  return(res)
}

check_occurrence_data <- function(species){
  filename = paste0("OccurrenceData/", species, ".csv")
  if(file.exists(filename)){
    res <- read.csv(filename, header = TRUE)
    return(res)
  } else {
    res <- get_occurrence_data(species)
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
  # Check the occurrence data, If there is an error there it catches it and writes an error in the output file
  tryCatch(res <- check_occurrence_data(species),
           error = function(errormessage){
             write.table(paste(c(species, rep("", 15), "ERROR while getting the occurrence data"), 
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

filter_n_closest_coordinate_ceiling <- function(n, occurrence_data, samplelocation){
  occurrence_data$distance <- pmax(abs(occurrence_data$Longitude - samplelocation$Longitude), 
                                   abs(occurrence_data$Latitude - samplelocation$Latitude))
  # Get the n-th smallest distance and round it up
  dist = ceiling(sort(occurrence_data$distance)[n])
  occurrence_data <- occurrence_data[occurrence_data$distance < dist,]
  return(occurrence_data)
}

sp_format <- function(coordinates){
  return(structure(as.numeric(c(coordinates["Longitude"], coordinates["Latitude"])), .Dim = 1:2))
}

find_shortest_route_in_sea <- function(samplelocation, occurrence_data, tr, row, filename){
  # Remove duplicates
  occurrence_data <- occurrence_data[!duplicated(occurrence_data),]
  # Remove samples taken further away than the closest point
  occurrence_data <- filter_on_distance(tr, samplelocation, occurrence_data)
  # Filter if there are more than 10 unique locations
  row$inrange <- nrow(occurrence_data)
  if(nrow(occurrence_data) > 10){
    occurrence_data <- filter_n_closest_coordinate_ceiling(10, occurrence_data, samplelocation)
  }
  # Save the number of points for which the distance will be calculated
  row$pointscalculated <- nrow(occurrence_data)
  # find the shortest route to every point through the sea
  paths <- sapply(1:nrow(occurrence_data), function(i) {
    path <- shortestPath(tr, sp_format(samplelocation), 
                         sp_format(occurrence_data[i,]), 
                         output = "SpatialLines")
    return(path)
  })
  # Find the closest location the point of sampling
  row$distance <- min(as.numeric(sapply(paths, function(x) geosphere::lengthLine(x))), na.rm=T)
  if(!file.exists(filename)){
    write.table(paste(c(names(row)), collapse = ","), file = filename, append = T, quote = F, sep = ",", col.names = F, row.names = F)
  }
  write.table(row, file = filename, append = T, quote = F, sep = ",", col.names = F, row.names = F)
}

plot_shortest_path <- function(path, SampleLocation, findLocations){
  ggplot() +
    geom_polygon(aes(x = long, y = lat, group = group), data = map_data("world")) +
    theme_void() +
    coord_cartesian(xlim = c(-25, 60), ylim = c(35,70)) +
    geom_point(aes(x= Longitude, y = Latitude), data = SampleLocation, col = "red") +
    geom_point(aes(x= Longitude, y = Latitude), data = findLocations, col = "black") +
    geom_path(aes(x=long,y=lat), color="blue", size = 2,
              data=fortify(SpatialLinesDataFrame(path, data = data.frame(ID = 1))))+
    geom_label_repel(aes(x= Longitude, y = Latitude, label = Observatory.ID), data = SampleLocation)
}

filter_on_distance <- function(tr, samplelocation, occurrence_data){
  distances <- distm(samplelocation, occurrence_data, fun = distVincentyEllipsoid)
  sea_dist <-  geosphere::lengthLine(shortestPath(tr, sp_format(samplelocation), 
                                                  sp_format(occurrence_data[which.min(distances),]), 
                                                  output = "SpatialLines"))
  return(occurrence_data <- occurrence_data[distances <= sea_dist,])
}