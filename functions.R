################################ Load Packages #################################
library("ggplot2")
library("leaflet")

################################## Functions ###################################

plotBar <- function(df) {
  plot <- ggplot(pivot_longer(df, names(df)[-1])) +
    theme(axis.text.x = element_text(angle = 15)) +
    geom_bar(aes(x = name, y = value, fill = name), stat = "identity") +
    scale_fill_manual("Sample Location", values = locCols)
  return(plot)
}

plotLeaflet <- function(df, locations, occurrence) {
  aicon = awesomeIcons(markerColor = locCols) 
  leaflet() %>%
    addProviderTiles(providers$Stamen.TonerLite,
    options = providerTileOptions(noWrap = TRUE)) %>%
      addAwesomeMarkers(data = locations[locations$Observatory.ID %in% names(df)[df[1,]>0],], ~Longitude, ~Latitude, label = ~Observatory.ID, icon = aicon) %>%
      addMarkers(data = occurrence, ~Longitude, ~Latitude, clusterOptions = markerClusterOptions())
}