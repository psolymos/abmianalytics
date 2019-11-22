library(shiny)
library(DT)
library(magick)

init <- "d:/tmp"
#init <- getwd()

ui <- fluidPage(

    titlePanel("WildTrax Picture Picker"),

    tabsetPanel(
        tabPanel("Settings",
            textInput("path", "Enter the path:", value=init, placeholder=init, width="100%"),
            verbatimTextOutput("pathvalue", placeholder = TRUE),
            checkboxGroupInput("filetypes", "Image types to consider:",
                choices=c("jpg", "png", "bmp", "tif"), selected=c("jpg", "png"), inline=TRUE),
            checkboxInput("casesens", "Should path and file names be case sensitive:"),
            hr(),
            dataTableOutput("dirindt")
        ),
        tabPanel("Pictures",
            #selectInput("batchsize", "Batch size:", c(16)),
            #uiOutput("pics")
            fluidRow(
                column(4, uiOutput("plot11")),
                column(4, uiOutput("plot12")),
                column(4, uiOutput("plot13"))
            ),
            fluidRow(
                column(4, uiOutput("plot21")),
                column(4, uiOutput("plot22")),
                column(4, uiOutput("plot23"))
            )
        ),
        tabPanel("Results",
            # dowload button
            # DT table with assessment
            tableOutput("table")
        )
    )
)


server <- function(input, output) {
    output$pathvalue <- renderText({
        input$path
    })
    get_filelist <- reactive({
        fl <- NULL
        for (i in input$filetypes)
            fl <- c(fl,
                list.files(input$path,
                    pattern = paste0("\\.", i, "$"),
                    ignore.case=!input$casesens))
        fl
    })
    output$dirindt <- renderDataTable({
        req(get_filelist())
        data.frame(get_filelist())
    })
    output$plot1 <- renderImage({
        req(get_filelist())
        fl <- get_filelist()
        list(src = file.path(input$path, fl[1]),
             alt = fl[1])
    })
    output$plot11 <- renderUI({
        i <- 1
        req(get_filelist())
        fl <- get_filelist()
        n <- length(fl)
        if (n < i)
            return(NULL)
        tagList(
            h4(fl[i]),
            renderImage({
                list(src = file.path(input$path, fl[i]),
                    alt = fl[1],
                    width="100%")
            }, deleteFile=FALSE)
        )
    })
    output$plot12 <- renderUI({
        i <- 2
        req(get_filelist())
        fl <- get_filelist()
        n <- length(fl)
        if (n < i)
            return(NULL)
        tagList(
            h4(fl[i]),
            renderImage({
                list(src = file.path(input$path, fl[i]),
                    alt = fl[1],
                    width="100%")
            }, deleteFile=FALSE)
        )
    })
    output$plot13 <- renderUI({
        i <- 3
        req(get_filelist())
        fl <- get_filelist()
        n <- length(fl)
        if (n < i)
            return(NULL)
        tagList(
            h4(fl[i]),
            renderImage({
                list(src = file.path(input$path, fl[i]),
                    alt = fl[1],
                    width="100%")
            }, deleteFile=FALSE)
        )
    })
    output$plot21 <- renderUI({
        i <- 4
        req(get_filelist())
        fl <- get_filelist()
        n <- length(fl)
        if (n < i)
            return(NULL)
        tagList(
            h4(fl[i]),
            renderImage({
                list(src = file.path(input$path, fl[i]),
                    alt = fl[1],
                    width="100%")
            }, deleteFile=FALSE)
        )
    })
    output$plot22 <- renderUI({
        i <- 5
        req(get_filelist())
        fl <- get_filelist()
        n <- length(fl)
        if (n < i)
            return(NULL)
        tagList(
            h4(fl[i]),
            renderImage({
                list(src = file.path(input$path, fl[i]),
                    alt = fl[1],
                    width="100%")
            }, deleteFile=FALSE)
        )
    })
    output$plot23 <- renderUI({
        i <- 6
        req(get_filelist())
        fl <- get_filelist()
        n <- length(fl)
        if (n < i)
            return(NULL)
        tagList(
            h4(fl[i]),
            renderImage({
                list(src = file.path(input$path, fl[i]),
                    alt = fl[1],
                    width="100%")
            }, deleteFile=FALSE)
        )
    })

}

shiny::shinyApp(ui, server)

#s:/reports/2018/images/birds/

