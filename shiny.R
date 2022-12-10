library("shiny")
library("ggplot2")
library("tidyr")
library("maps")
library("rlist")
library("rgbif")
library("leaflet")
# library("rentrez")
library("RColorBrewer")

#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("functions.R")
df <- readRDS("inputs/BOLDigger_Species_Location.rds")
locations <- read.csv("inputs/Coordinates.csv")

# locations$Observatory.ID
locCols <- brewer.pal(nrow(locations), "Set3")
names(locCols) <- locations$Observatory.ID

ui <- pageWithSidebar(
  # App title ----
  titlePanel("GBIF occurrence"),
  # Sidebar panel for inputs ----
  sidebarPanel(
    selectInput("species", "Species", df$Specieslist),
    # radioButtons("line", "Show shortest path through sea", c("No" = FALSE, "Yes" = TRUE)),
    # radioButtons("visualization", "Show occurrences as", c("Points" = geom_point, "Hexagons" = geom_hex),
    # radioButtons("allLoc", "Show all sample locations",c("Yes" = TRUE, "No" = FALSE)),
    width = 2,
    tableOutput("sptable")
  ),
  # Main panel for displaying outputs ----
  mainPanel(
    plotOutput("barplot"),
    leafletOutput("mymap"),
    width = 10
  )
)
server <- function(input, output){
  output$species <- renderText(input$species)
  dfSp <- eventReactive(input$species, df[df$Specieslist == input$species, ])
  infoSp <- eventReactive(input$species, t(name_backbone(input$species)[,c("usageKey", "scientificName", "kingdom", "phylum", "order", "family", "class")]))
  # NCBI <- eventReactive(input$species, entrez_search(db="taxonomy", term=paste0("(", input$species, "[ORGN]) AND Species[RANK]"))$ids)
  # infoSp <- eventReactive(input$species, rbind(infoSp, data.frame("NCBI ID" = NCBI)))
  occurrence <- eventReactive(input$species, {
    read.csv(paste0("OccurrenceData/", input$species, ".csv"))
  })
  output$barplot <- renderPlot(plotBar(dfSp()))
  output$sptable <- renderTable(infoSp(), rownames = TRUE, colnames = FALSE)
  output$mymap <- renderLeaflet(plotLeaflet(dfSp(), locations, occurrence()))
}


options(shiny.host = "0.0.0.0")
options(shiny.port = 80)
shinyApp(ui, server)


# species = "Homo sapiens"
# infoSp <- t(name_backbone(species)[,c("usageKey", "scientificName", "kingdom", "phylum", "order", "family", "class")])
# NCBI <- eventReactive(input$species, entrez_search(db="taxonomy", term=paste0("(", input$species, "[ORGN]) AND Species[RANK]"))$ids
# infoSp <- rbind(infoSp, data.frame("NCBI ID" = NCBI))