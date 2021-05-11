metaboVisUI <- function(id){
    ns <- NS(id)
    fluidPage(
        fluidRow(
            column(width = 4,
                selectInput(inputId = ns("object"), label = "Metabolomics object:", choices = ""),
                selectInput(inputId = ns("files"), label = "Files:", choices = "", multiple = TRUE),
                selectInput(inputId = ns("chromtype"), label = "Select chromatogrma type:", choices = c("max", "sum"))
            ),
            column(width = 4,
                selectInput(inputId = ns("groupVar"), label = "Grouping variable to make comparisons:",
                                        choices = ""),
                numericInput(inputId = ns("mz"), label = "mz value", value = 500, min = 0, max = 3000, step = 0.001),
                numericInput(inputId = ns("ppm"), label = "ppm error", value = 30, min = 1, max = 500)
            ),
        ),
        fluidRow(
            column(width = 8,
                sliderInput(inputId = ns("rt"), label = "RT range", min = 0,
                            max = 1000, value = c(500,750), step = 1)
            )
        ),
        fluidRow(
            plotlyOutput(ns("chromplot"))
        )
    )
}
