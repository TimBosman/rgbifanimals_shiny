################################ Load Packages #################################
library("ggplot2")
library("shiny")
library("tidyr")
library("maps")
library("rlist")
library("leaflet")
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
  # load all the data sets if the species changes
  df_sp <- eventReactive(input$species, df[df$Specieslist == input$species, ])
  info_sp <- eventReactive(input$species, data_gbif(input$species))
  occurrence <- eventReactive(input$species, data_occurence(input$species))

  # define all the tables and plots here
  output$barplot <- renderPlot(plot_bar(df_sp(), loc_cols))
  output$sptable <- renderTable(info_sp(), rownames = TRUE, colnames = FALSE)
  output$mymap <- renderLeaflet(plot_leaflet(df_sp(), locations, occurrence(),
    loc_cols))
}

###################### Run the Shiny app on localhost:80 #######################
options(shiny.host = "0.0.0.0")
options(shiny.port = 80)
shinyApp(ui, server)
