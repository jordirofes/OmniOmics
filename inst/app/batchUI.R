batchUI <- function(id){
    ns <- NS(id)
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
                                choices = "")
                ),
                column(width = 4,
                    radioButtons(inputId = ns("phenoVar2"), label = "Second covariate (optional):",
                                choices = "")
                )
            , ns = ns)
        )
    )
}
