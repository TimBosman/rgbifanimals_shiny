################################ Load Packages #################################
library("ggplot2")
library("shiny")
library("tidyr")
library("maps")
library("rlist")
library("rgbif")
library("leaflet")
# library("rentrez")
library("RColorBrewer")

######################### Load functions and datasets ##########################
source("functions.R")
df <- readRDS("inputs/BOLDigger_Species_Location.rds")
locations <- read.csv("inputs/Coordinates.csv")

loc_cols <- brewer.pal(nrow(locations), "Set3")
names(loc_cols) <- locations$Observatory.ID

################################ Define the GUI ################################
ui <- pageWithSidebar(
  # App title ----
  titlePanel("GBIF occurrence"),
  # Sidebar panel for inputs ----
  sidebarPanel(
    selectInput("species", "Species", df$Specieslist),
    # radioButtons("line", "Show shortest path through sea",
      # c("No" = FALSE, "Yes" = TRUE)),
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

############################## Define the backend ##############################
server <- function(input, output) {
  output$species <- renderText(input$species)
  df_sp <- eventReactive(input$species, df[df$Specieslist == input$species, ])
  info_cols <- c("usageKey", "scientificName", "kingdom", "phylum", "order",
    "family", "class")
  info_sp <- eventReactive(input$species,
    t(name_backbone(input$species)[, info_cols]))
  # NCBI <- eventReactive(input$species, entrez_search(db="taxonomy",
  #   term=paste0("(", input$species, "[ORGN]) AND Species[RANK]"))$ids)
  # info_sp <- eventReactive(input$species,
  #   rbind(infoSp, data.frame("NCBI ID" = NCBI)))
  occurrence <- eventReactive(input$species, {
    read.csv(paste0("OccurrenceData/", input$species, ".csv"))
  })
  output$barplot <- renderPlot(plot_bar(df_sp(), loc_cols))
  output$sptable <- renderTable(info_sp(), rownames = TRUE, colnames = FALSE)
  output$mymap <- renderLeaflet(plot_leaflet(df_sp(), locations, occurrence(),
    loc_cols))
}

###################### Run the Shiny app on localhost:80 #######################
options(shiny.host = "0.0.0.0")
options(shiny.port = 80)
shinyApp(ui, server)
