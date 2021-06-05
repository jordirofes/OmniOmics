batchUI <- function(id){
    ns <- NS(id)
    tabsetPanel(
        tabPanel(title = "Batch Correction",
            fluidPage(
                fluidRow(
                    column(width = 4,
                           selectInput(inputId = ns("object"), label = "Select an object to process:",
                                       choices = "")
                    ),
                    column(width = 4,
                           selectInput(inputId = ns("batchNorm"), label = "Select a normalization method:",
                                       choices = c("Covariate", "Injection Order"))
                    ),
                    column(width = 4,
                           actionButton(inputId = ns("batchBut"), label = "Begin batch correction:")
                    )
                ),
                fluidRow(
                    conditionalPanel("input.batchNorm == 'Injection Order'",
                        column(width = 4,
                                radioButtons(inputId = ns("groupVar"), label = "QC/Sample/Blank variable:",
                                             choices = "")
                        ),
                        column(width = 4,
                                radioButtons(inputId = ns("ordVar"), label = "Order Variable:",
                                             choices = "")
                        ),
                        column(width = 4,
                                radioButtons(inputId = ns("batchVar"), label = "Batch Variable:",
                                             choices = "")
                        ),
                        column(width = 4,
                                selectInput(inputId = ns("qcname"), label = "Select QC name", choices = "")
                        )
                    , ns = ns),
                    conditionalPanel("input.batchNorm == 'Covariate'",
                        column(width = 4,
                                radioButtons(inputId = ns("phenoVar"), label = "Covariate:",
                                             choices = "", selected = character(0))
                        ),
                        column(width = 4,
                                radioButtons(inputId = ns("phenoVar2"), label = "Second covariate (optional):",
                                             choices = "",selected = character(0))
                        )
                    , ns = ns)
                )
            )
        ),
        tabPanel(title = "Batch Assessment",
            fluidPage(
                fluidRow(
                    column(width = 4,
                           selectInput(inputId = ns("object2"), label = "Select an object to process:",
                                       choices = "")
                    ),
                    column(width = 4,
                           selectInput(inputId = ns("batchNorm2"), label = "Select a normalization method:",
                                       choices = c("Covariate", "Injection Order"))
                    )
                ),
                conditionalPanel("input.batchNorm2 == 'Covariate'",
                    fluidRow(
                        column(width = 4,
                            selectInput(inputId = ns("phenoVar3"), label = "Covariate:",
                                        choices = "", multiple = TRUE)
                        ),
                        column(width = 4,
                            numericInput(inputId = ns("thr"), label = "% minimum variability by selected principal components:",
                                        value = 0.3, min = 0, max = 1, step = 0.01)
                        )
                    ),
                    fluidRow(
                        column(width = 12,
                            plotOutput(ns("covplot"))
                        )
                    )
                , ns = ns),
                conditionalPanel("input.batchNorm2 == 'Injection Order'",
                    fluidRow(
                        column(width = 4,
                            selectInput(inputId = ns("groupVar2"), label = "QC/Sample/Blank variable:",
                                        choices = "")
                        )
                    ),
                    fluidRow(
                        column(width = 12,
                            plotlyOutput(ns("injplot"))
                        )
                    )
                , ns = ns)
            )
        )
    )
}
