objectDataUI <- function(id, objectList = ""){
    ns <- NS(id)
    fluidPage(
        fluidRow(
            column(width = 4,
                    selectInput(inputId = ns("object"), label = "Select an object:",
                                choices = "")
            )
        ),
        fluidRow(
            column(width = 12, offset = 0,
                    verbatimTextOutput(outputId = ns("dt_info"), placeholder = TRUE)
            )
        ),
        fluidRow(
            column(width = 12,
                    DTOutput(outputId = ns("extra_dt"))
            )
        )
    )
}
