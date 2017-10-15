source("globals.R")

ui <- fluidPage(
  titlePanel("Human Footprint in Alberta"),
    wellPanel(fluidRow(
      column(4,
        checkboxGroupInput("checkRegions", label = "Natural Regions",
          choices = list("Grassland"=1, "Parkland"=2, "Foothills"=3, "Boreal"=4,
          "Canadian Shield"=5, "Rocky Mountain"=6),
          selected = 1:6)
      ),
      column(4,
        checkboxGroupInput("checkSectors", label = "Industrial Sectors",
          choices = list("Agriculture"=1, "Forestry"=2, "Energy"=3,
          "Urban"=4, "Transportation"=5, "Other"=6),
          selected = 1:6)
      ),
      column(4,
        radioButtons("byReg", label = "Summarize by",
          choices = list("Industrial Sector" = 1, "Natural Region" = 2),
          selected = 1)
      )
    )),
    fluidRow(
      column(6,
        htmlOutput("myChart"),
        p("Percent human footprint according to selection. Based on yearly veryfied human footprint in 3 x 7 km rectangles at ABMI site locations. Hover over the lines to see values, use Edit button to change chart settings.")
      ),
      column(6,
#        plotOutput("myMap"),
        leafletOutput("myMap"),
        p("The map shows percent human footprint according to selection. Based on 2014 human footprint inventory at 10 x 10 km scale.")
      )
    )
)

server <- function(input, output) {
    output$myChart <- renderGvis({
        get_gplot(as.integer(input$checkRegions),
            as.integer(input$checkSectors), input$byReg == 2)
    })
#    output$myPlot <- renderPlot({
#        get_plot(as.integer(input$checkRegions),
#            as.integer(input$checkSectors), input$byReg == 2)
#    })
#    output$myMap <- renderPlot({
#        get_rmap(as.integer(input$checkRegions), as.integer(input$checkSectors))
#    })
    output$myMap <- renderLeaflet({
        get_mape(as.integer(input$checkRegions), as.integer(input$checkSectors))
    })
}

shiny::shinyApp(ui, server)

# regions/sectors < 2
