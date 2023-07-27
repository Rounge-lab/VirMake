require("shiny")
require("shinydashboard")
require("shinythemes")

UI <- fluidPage(theme = shinytheme("cerulean"),
    navbarPage(
        "VirMake results",
        tabPanel("QC"),
        tabPanel("Assembly"),
        tabPanel("Identification"),
        tabPanel("Mapping"),
        tabPanel("Taxonomy"),
        tabPanel("Functional analysis"),
    )

)

server <- function(input, output) {
    output$txtout <- renderText({
        paste(input$txt1, input$txt2, sep = " ")
    })
}

shinyApp(ui = UI, server = server)
