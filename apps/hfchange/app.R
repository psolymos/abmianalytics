source("globals.R")

ui <- fluidPage(
  titlePanel("Human Footprint in Alberta"),
    fluidRow(
      column(2,
        wellPanel(
            h3("Settings"),
            radioButtons("byReg", label = "Summarize by",
                choices = list("Industrial Sector" = 1, "Natural Region" = 2),
                selected = 1),
            checkboxGroupInput("checkRegions", label = "Natural Regions",
                choices = list("Grassland"=1, "Parkland"=2, "Foothills"=3, "Boreal"=4,
                "Canadian Shield"=5, "Rocky Mountain"=6),
                selected = 1:6),
            checkboxGroupInput("checkSectors", label = "Natural Regions",
                choices = list("Agriculture"=1, "Forestry"=2, "Energy"=3,
                "Urban"=4, "Transportation"=5, "Other"=6),
                selected = 1:6)
        )
      ),
      column(4, plotOutput("myMap")),
      column(6, htmlOutput("myChart"))
    )
)

server <- function(input, output) {
    output$myChart <- renderGvis({
        get_gplot(as.integer(input$checkRegions),
            as.integer(input$checkSectors), input$byReg == 2)
    })
    output$myPlot <- renderPlot({
        get_plot(as.integer(input$checkRegions),
            as.integer(input$checkSectors), input$byReg == 2)
    })
    output$myMap <- renderPlot({
        get_map(input$checkRegions)
    })
}

shiny::shinyApp(ui, server)

# regions/sectors < 2
