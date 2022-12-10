################################ Load Packages #################################
library("ggplot2")
library("leaflet")
library("tidyr")

################################## Functions ###################################

plot_bar <- function(df, loc_cols) {
  plot <- ggplot(pivot_longer(df, names(df)[-1])) +
    theme(axis.text.x = element_text(angle = 15)) +
    geom_bar(aes(x = name, y = value, fill = name), stat = "identity") +
    scale_fill_manual("Sample Location", values = loc_cols)
  return(plot)
}

plot_leaflet <- function(df, locations, occurrence, loc_cols) {
  locations <- locations[locations$Observatory.ID %in% names(df)[df[1, ] > 0], ]
  aicon <- awesomeIcons(markerColor = loc_cols)
  leaflet() %>%
    addProviderTiles(providers$Stamen.TonerLite,
    options = providerTileOptions(noWrap = TRUE)) %>%
      addAwesomeMarkers(data = locations, ~Longitude, ~Latitude,
        label = ~Observatory.ID, icon = aicon) %>%
      addMarkers(data = occurrence, ~Longitude, ~Latitude,
        clusterOptions = markerClusterOptions())
}