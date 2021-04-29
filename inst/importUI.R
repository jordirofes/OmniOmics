importUI <- function(id){
    ns <- NS(id)
    fluidPage(
        fluidRow(
            column(width = 4,
                    switchInput(inputId = ns("header"),onLabel = "TRUE",
                                offLabel = "FALSE", value = FALSE, label = "Header"),
                    radioButtons(inputId = ns("sep"), label = "Separator:",
                                choiceNames = c("Comma", "Semicolon", "Tab", "Empty Space"),
                                choiceValues = c(",", ";", "\t", " "))
            ),
            column(width = 4, offset = 0,
                    fileInput(inputId = ns("pd"), label = "Load your Files:",
                            buttonLabel = "Load", multiple = TRUE)

            )
        ),
        fluidRow(
            column(width = 12,
                DT::DTOutput(ns("phenoTable"))
            )
        ),
        uiOutput(ns("dt_upload_ui"))
    )
}
