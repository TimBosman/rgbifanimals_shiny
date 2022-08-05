################################### Load Packages ############################################
library("ggplot2")

################################### Functions ################################################

plotMap <- function(df, locations, occurrence, size = ""){
  plot <- ggplot()+
    geom_polygon(aes(x = long, y = lat, group = group), data = map_data("world")) +
    theme_void() + 
    geom_hex(data = occurrence, aes(x = Longitude, y = Latitude), alpha = 0.8, col = "White") +
    scale_fill_gradient(low = "#8CDA8F", high = "#0D1F0B") +
    geom_point(data = locations, aes(x = Longitude, y = Latitude, col = Observatory.ID), pch = 16, size = 5)+
    scale_color_manual("Sample Location", values = locCols) 
    if(size == "europe"){
      plot <- plot + coord_cartesian(xlim = c(-25, 60), ylim = c(35,70))
    }
    if(size == "eastsea"){
      plot <- plot + coord_cartesian(xlim = c(0, 30), ylim = c(50,65))
  }
  print(plot)
}