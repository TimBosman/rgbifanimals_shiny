library("shiny")
library("ggplot2")
library("tidyr")
library("maps")
library("rlist")
library("rgbif")
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
    radioButtons("mapsize", "Map size", c("World" = "world", "Europe" = "europe", "East Sea" = "eastsea"), selected = "europe"),
    # radioButtons("line", "Show shortest path through sea", c("No" = FALSE, "Yes" = TRUE)),
    # radioButtons("visualization", "Show occurrences as", c("Points" = geom_point, "Hexagons" = geom_hex),
    # radioButtons("allLoc", "Show all sample locations",c("Yes" = TRUE, "No" = FALSE)),
    width = 2,
    tableOutput("sptable")
  ),
  # Main panel for displaying outputs ----
  mainPanel(
    plotOutput("barplot"),
    plotOutput("mapplot"),
    width = 10
  )
)
server <- function(input, output){
  output$species <- renderText(input$species)
  dfSp <- eventReactive(input$species, df[df$Specieslist == input$species, ])
  infoSp <- eventReactive(input$species, t(name_backbone(input$species)[,c("usageKey", "scientificName", "kingdom", "phylum", "order", "family", "class")]))
  occurrence <- eventReactive(input$species, {
    read.csv(paste0("OccurrenceData/", input$species, ".csv"))
  })
  output$barplot <- renderPlot(plotBar(dfSp()))
  output$mapplot <- renderPlot(plotMap(dfSp(), locations, occurrence(), input$mapsize))
  output$sptable <- renderTable(infoSp(), rownames = TRUE, colnames = FALSE)
}

shinyApp(ui,server, options = c(port = 80))
