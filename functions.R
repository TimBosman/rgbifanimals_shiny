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
    scale_fill_manual("Sample Location", values = loc_cols) +
    xlab("Sample location") +
    ylab("Number of sequences found")
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
      icon = sample_icons, group = "Sample locations",
      label = ~Observatory.ID) %>%
    addAwesomeMarkers(data = occurrence, ~Longitude, ~Latitude,
      icon = occurrence_icons, group = "GBIF locations",
      clusterOptions = markerClusterOptions()) %>%
    addLayersControl(overlayGroups = c("Sample locations", "GBIF locations"),
    options = layersControlOptions(collapsed = FALSE)
  )
}

################################ Data functions ################################
data_gbif <- function(species) {
  info_cols <- c("usageKey", "kingdom", "phylum", "class", "order",
    "family", "scientificName")
  # Retrieve GBIF information
  species_info <- t(name_backbone(species)[, info_cols])
  # Retrieve NCBI information
  NCBIquery <- paste0("(", species, "[ORGN]) AND Species[RANK]")
  NCBI <- entrez_search(db = "taxonomy", term = NCBIquery)$ids
  common_name <- entrez_summary(db="taxonomy", id=NCBI)$commonname
  # Combine the information
  species_info <- rbind(species_info, "NCBI ID" = NCBI)
  species_info <- rbind("Common name" = common_name, species_info)
  return(species_info)
}

data_occurence <- function(species) {
  return(read.csv(paste0("OccurrenceData/", species, ".csv")))
}