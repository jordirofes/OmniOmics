importUI <- function(id){
    ns <- NS(id)
    fluidPage(
        tabsetPanel(
            tabPanel(title = "Raw Data",
                fluidRow(
                column(width = 4,
                    switchInput(inputId = ns("header"),onLabel = "TRUE",
                                offLabel = "FALSE", value = FALSE, label = "Header"),
                    radioButtons(inputId = ns("sep"), label = "Separator:",
                                choiceNames = c("Comma", "Semicolon", "Tab", "Empty Space"),
                                choiceValues = c(",", ";", "\t", " "))
                ),
                column(width = 4, offset = 0,
                    fileInput(inputId = ns("pd"), label = "Load your PhenoData:",
                            buttonLabel = "Load", multiple = TRUE)

                )
                ),
                fluidRow(
                    column(width = 12,
                        DT::DTOutput(ns("phenoTable"))
                    )
                ),
                uiOutput(ns("dt_upload_ui"))
            ),
            tabPanel(title = "Files",
                fluidRow(
                    column(width = 4,
                        shinyFilesButton(ns("dt_path2"),style = "margin-top: 25px;", label="File select", title="Please select a file", multiple=TRUE),
                        actionButton(ns("load_data2"), label = "Load", style = "margin-top: 25px;"),
                        verbatimTextOutput(outputId = ns("path_names2"), placeholder = TRUE)
                    ),
                    column(width = 4,
                        textInput(inputId = ns("fileName2"), value = "experiment1", label = "Enter a name for the loaded object:")
                    )
                ),
                fluidRow(
                    column(width = 4,
                        shinySaveButton(ns("folder"),style = "margin-top: 25px;", label="Save data:", title="Select a folder to save the data:"),
                        actionButton(ns("save_data"), label = "Save", style = "margin-top: 25px;"),
                        verbatimTextOutput(outputId = ns("dir_name"), placeholder = TRUE)
                    ),
                    column(width = 4,
                        selectInput(inputId = ns("object"), label = "Select objects to save:",
                                    choices = "", multiple = TRUE)
                    )
                )
            )
        )
    )
}
