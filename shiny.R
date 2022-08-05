library("shiny")
library("ggplot2")
library("tidyr")
library("maps")
library("rlist")
library("RColorBrewer")

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
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
    selectInput("species", "Species", df$Specieslist, width = "100pt"),
    radioButtons("mapsize", "Map size", c("World" = "world", "Europe" = "europe", "East Sea" = "eastsea"), width = "100pt", selected = "europe"),
    # radioButtons("line", "Show shortest path through sea", c("No" = FALSE, "Yes" = TRUE), width = "100pt"),
    # radioButtons("visualization", "Show occurrences as", c("Points" = geom_point, "Hexagons" = geom_hex), width = "100pt"),
    # radioButtons("allLoc", "Show all sample locations",c("Yes" = TRUE, "No" = FALSE))
  ),
  # Main panel for displaying outputs ----
  mainPanel(
    plotOutput("barplot", width = "380pt"),
    plotOutput("mapplot", width = "380pt")
  )
)
server <- function(input, output){
  output$species <- renderText(input$species)
  dfSp <- eventReactive(input$species, {
    df[df$Specieslist == input$species, ]
  }) 
  occurrence <- eventReactive(input$species, {
    read.csv(paste0("OccurrenceData/", input$species, ".csv"))
  })
  output$barplot <- renderPlot(ggplot(pivot_longer(dfSp(), names(dfSp())[-1])) + 
                                 theme(axis.text.x = element_text(angle = 15)) +
                                 geom_bar(aes(x = name, y = value, fill = name), stat = "identity") +
                                 scale_fill_manual("Sample Location", values = locCols))
  output$mapplot <- renderPlot(plotMap(dfSp(), locations, occurrence(), input$mapsize)
  )
}

shinyApp(ui,server, options = c(port = 80))