################################### Load Packages ############################################
library("ggplot2")

################################### Functions ################################################

plotMap <- function(df, locations, occurrence, size = ""){
  plot <- ggplot()+
    geom_polygon(aes(x = long, y = lat, group = group), data = map_data("world")) +
    theme_void() + 
    geom_hex(data = occurrence, aes(x = Longitude, y = Latitude), alpha = 0.8, col = "White") +
    scale_fill_gradient(low = "#8CDA8F", high = "#0D1F0B") +
    geom_point(data = locations[locations$Observatory.ID %in% names(df)[df[1,]>0],], 
               aes(x = Longitude, y = Latitude, col = Observatory.ID), pch = 16, size = 5)+
    scale_color_manual("Sample Location", values = locCols) 
    if(size == "europe"){
      plot <- plot + coord_cartesian(xlim = c(-25, 60), ylim = c(35,70))
    }
    if(size == "eastsea"){
      plot <- plot + coord_cartesian(xlim = c(0, 30), ylim = c(50,65))
  }
  return(plot)
}

plotBar<- function(df){
  plot <- ggplot(pivot_longer(df, names(df)[-1])) + 
    theme(axis.text.x = element_text(angle = 15)) +
    geom_bar(aes(x = name, y = value, fill = name), stat = "identity") +
    scale_fill_manual("Sample Location", values = locCols)
  return(plot)
}