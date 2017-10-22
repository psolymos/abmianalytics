source("globals.R")

ui <- fluidPage(
  titlePanel("Human Footprint in Alberta"),
    fluidRow(
      column(6,
        htmlOutput("myChart"),
        p("Percent human footprint according to selection. Based on yearly veryfied human footprint in 3 x 7 km rectangles at ABMI site locations. Hover over the lines to see values, use Edit button to change chart settings.")
      ),
      column(6,
        leafletOutput("myMap", width="100%", height="450"),
        p("The map shows percent human footprint according to selection. Based on 2014 human footprint inventory at 4 x 4 km scale.")
      )
    ),
    wellPanel(fluidRow(
      column(3,
        radioButtons("byReg", label = "Summarize by",
          choices = list("Industrial Sector" = 1, "Region" = 2),
          selected = 1)
      ),
      column(3,
        checkboxGroupInput("checkSectors", label = "Industrial Sectors",
          choices = list(
            "Agriculture"=1,
            "Forestry"=2,
            "Energy"=3,
            "Urban"=4,
            "Transportation"=5,
            "Other"=6),
            selected = 1:6)
      ),
      column(3,
        radioButtons("whichReg", label = "Regions",
          choices = list(
            "Natural Regions" = 1,
            "Natural Subregions" = 2,
            "Land Use Framework Regions" = 3),
          selected = 1)
      ),
      column(3,
        conditionalPanel(condition="input.whichReg == '1'",
          checkboxGroupInput("checkRegionsNR", label = "Natural Regions",
            choices = list(
                "Boreal"=1,
                "Canadian Shield"=2,
                "Foothills"=3,
                "Grassland"=4,
                "Parkland"=5,
                "Rocky Mountain"=6),
            selected = 1:6)),
        conditionalPanel(condition="input.whichReg == '2'",
          checkboxGroupInput("checkRegionsNSR", label = "Natural Subregions",
            choices = list(
                "Alpine"=1,
                "Athabasca Plain"=2,
                "Boreal Subarctic"=3,
                "Central Mixedwood"=4,
                "Central Parkland"=5,
                "Dry Mixedgrass"=6,
                "Dry Mixedwood"=7,
                "Foothills Fescue"=8,
                "Foothills Parkland"=9,
                "Kazan Uplands"=10,
                "Lower Boreal Highlands"=11,
                "Lower Foothills"=12,
                "Mixedgrass"=13,
                "Montane"=14,
                "Northern Fescue"=15,
                "Northern Mixedwood"=16,
                "Peace-Athabasca Delta"=17,
                "Peace River Parkland"=18,
                "Subalpine"=19,
                "Upper Boreal Highlands"=20,
                "Upper Foothills"=21),
            selected = 1:21)),
        conditionalPanel(condition="input.whichReg == '3'",
          checkboxGroupInput("checkRegionsLUF", label = "Land Use Framework Regions",
            choices = list(
                "Lower Athabasca"=1,
                "Lower Peace"=2,
                "North Saskatchewan"=3,
                "Red Deer"=4,
                "South Saskatchewan"=5,
                "Upper Athabasca"=6,
                "Upper Peace"=7),
            selected = 1:7))
      )
    ))
)

server <- function(input, output) {
    output$myChart <- renderGvis({
        rval <- switch(as.character(input$whichReg),
            "1"=as.integer(input$checkRegionsNR),
            "2"=as.integer(input$checkRegionsNSR),
            "3"=as.integer(input$checkRegionsLUF))
        get_gplot(
            r=rval,
            c=as.integer(input$checkSectors),
            byregion=input$byReg == 2,
            which_region=as.integer(input$whichReg))
    })
    output$myMap <- renderLeaflet({
        rval <- switch(as.character(input$whichReg),
            "1"=as.integer(input$checkRegionsNR),
            "2"=as.integer(input$checkRegionsNSR),
            "3"=as.integer(input$checkRegionsLUF))
        get_mape(
            r=rval,
            c=as.integer(input$checkSectors),
            which_region=as.integer(input$whichReg))
    })
}

shiny::shinyApp(ui, server)

