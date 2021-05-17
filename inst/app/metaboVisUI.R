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
                switchInput(inputId = ns("byrange"), label = "mz by range")
            ),
            column(width = 4,
                switchInput(inputId = ns("log"), label = "Logscaled")
            )
        ),
        fluidRow(
            column(width = 6,
                sliderInput(inputId = ns("rt"), label = "RT range", min = 0,
                            max = 1000, value = c(500,750), step = 1)
            ),
            conditionalPanel("input.byrange",
                column(width = 6,
                    sliderInput(inputId = ns("mz1"), label = "MZ range", min = 0,
                                max = 1000, value = c(200,300), step = 0.001)
                )
            , ns = ns),
            conditionalPanel("input.byrange == false",
                column(width = 6,
                    numericInput(inputId = ns("mz2"), label = "mz value", value = 500, min = 0, max = 3000, step = 0.001),
                    numericInput(inputId = ns("ppm"), label = "ppm error", value = 30, min = 1, max = 500)
                )
            , ns = ns)
        ),
        fluidRow(
            plotlyOutput(ns("chromplot"))
        )
    )
}
