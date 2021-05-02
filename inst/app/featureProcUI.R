featureProcUI <- function(id){
    ns <- NS(id)
    fluidPage(
        fluidRow(
            column(width = 4,
                selectInput(inputId = ns("object"), label = "Select an object to process:",
                            choices = "")
            ),
            column(width = 4,
                selectInput(inputId = ns("omic"), label = "Select an omic:",
                            choices = c("Metabolomics", "Transcriptomics"))
            ),
        ),
        conditionalPanel("input$omic == 'Metabolomics'",
            fluidRow(
                column(width = 4,
                    switchInput(ns("blankFilt"), label = "Blank Filter"),
                    switchInput(ns("cvFilt"), label = "Filter by variation coefficient"),
                    switchInput(ns("naFilt"), label = "Filter by NA ratio"),
                ),
                column(width = 4,
                    switchInput(ns("metvarFilt"), label = "Variance Filter"),
                    switchInput(ns("intFilt"), label = "Intensity Filter"),
                    switchInput(ns("sampleFilt"), label = "Sample Filter")
                ),
                column(width = 4,
                    switchInput(ns("ism0Filt"), label = "Filter non M0 features"),
                    switchInput(ns("annoFilt"), label = "Filter non annotated features")
                ),
            ),
            fluidRow(
                column(width = 4,
                    radioButtons(inputId = ns("groupVar"), label = "QC/Sample/Blank variable:",
                                choices = ""),
                    selectInput(inputId = ns("blankname"), label = "Select a blank name", choices = ""),
                    selectInput(inputId = ns("qcname"), label = "Select the QC name", choices = ""),
                    selectInput(inputId = ns("filtername"), label = "Select the sample name to filter", choices = "")
                ),
                column(width = 4,
                    conditionalPanel("input$blankFilt",
                        numericInput(inputId = ns("blankRatio"), label = "Select a fold-change between sample and blank:",
                                    value = 1, min = 0, max = 100, step = 0.1)

                    , ns = ns),
                    conditionalPanel("input$cvFilt",
                        numericInput(inputId = ns("maxRDS"),
                                    label = "Select a maximum QC RDS % threshold::",
                                    value = 20, min = 0, max = 100)
                    , ns = ns),
                    conditionalPanel("input$naFilt",
                        numericInput(inputId = ns("nathr"), label = "NA ratio threhshold",
                                    value = 0.4, min = 0, max = 1, step = 0.01),
                        selectInput(inputId = ns("naMethod"), label = "Select a NA ratio method:",
                                    choices = c(WithinQC = "QC", WithinSampleClass = "within",
                                                AcrossAllSamples = "across"))
                    , ns = ns)
                ),
                column(width = 4,
                    conditionalPanel("input$metvarFilt",
                        selectInput(inputId = ns("varfun"), label = "Select a variance function:",
                                    choices = c("IQR", "CV", "SD")),
                        switchInput(inputId = ns())
                    ),
                    conditionalPanel("input$intFilt",

                    ),
                    conditionalPanel("input$sampleFilt",

                    ),
                )
            )
        , ns = ns),
        conditionalPanel("input$omic == 'Transcriptomics'",
            fluidRow(
                column(width = 4,
                    switchInput(inputId = ns("entrez"), label = "Remove non entrez genes and duplicated")
                ),
                column(width = 4,
                    switchInput(inputId = ns("varfilt"), label = "Variance Filter")
                ),
                column(width = 4,
                    switchInput(inputId = ns("compFilter"), label = "Remove non-significant features")
                )
            ),
            fluidRow(
                conditionalPanel("input$varfilt",
                    column(width = 4,
                        selectInput(inputId = ns("varfun"), label = "Select the variance function to apply:"),
                        numericInput(inputId = ns("quantvar"), label = "Select the quantile to eliminate:",
                                    value = 0.6, min = 0, max = 1, step = 0.01)
                    )
                ,ns = ns),
                conditionalPanel("input$compFilter",
                    column(width = 4,
                        selectInput(inputId = ns("compTable"),
                                    label = "Select the comparison table to use for filtering:",
                                    choices = ""),
                        switchInput(inputId = ns("padj"), label = "Adjusted p-value?"),
                        numericInput(inputId = ns("pthr"), label = "P-value threshold",
                                    value = 0.05, min = 0, max = 1, step = 0.01),
                        numericInput(inputId = ns("fcthr"), label = "log2(Fold-Change) threshold",
                                     value = 1, min = 0, max = 100, step = 1)
                    )
                )
            )
        , ns = ns)
    )
}
