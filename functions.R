################################ Load Packages #################################
library("ggplot2")
library("leaflet")
library("tidyr")
library("rentrez")
library("rgbif")

############################### Plot  functions ################################

plot_bar <- function(df, loc_cols) {
  plot <- ggplot(pivot_longer(df, names(df)[-1])) +
    theme(axis.text.x = element_text(angle = 15)) +
    geom_bar(aes(x = name, y = value, fill = name), stat = "identity") +
    scale_fill_manual("Sample Location", values = loc_cols)
  return(plot)
}

plot_leaflet <- function(df, locations, occurrence, loc_cols) {
  locations <- locations[locations$Observatory.ID %in% names(df)[df[1, ] > 0], ]
  occurrence_icons <- awesomeIcons(icon = "bookmark", markerColor = "black")
  sample_icons <- awesomeIcons(icon = "flag", markerColor = "green")
  leaflet() %>%
    addProviderTiles(providers$Stamen.TonerLite,
    options = providerTileOptions(noWrap = TRUE)) %>%
      addAwesomeMarkers(data = locations, ~Longitude, ~Latitude,
        label = ~Observatory.ID, icon = sample_icons) %>%
      addAwesomeMarkers(data = occurrence, ~Longitude, ~Latitude,
        clusterOptions = markerClusterOptions(), icon = occurrence_icons)
}

################################ Data functions ################################
data_gbif <- function(species) {
  info_cols <- c("usageKey", "scientificName", "kingdom", "phylum", "order",
    "family", "class")
  species_info <- t(name_backbone(species)[, info_cols])

  # NCBI <- eventReactive(species, entrez_search(db="taxonomy",
  #   term=paste0("(", species, "[ORGN]) AND Species[RANK]"))$ids)
  # info_sp <- eventReactive(species,
  #   rbind(infoSp, data.frame("NCBI ID" = NCBI)))

  return(species_info)
}

data_occurence <- function(species) {
  return(read.csv(paste0("OccurrenceData/", species, ".csv")))
}