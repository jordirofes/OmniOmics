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
                            choices = c("Metabolomics", "Transcriptomics", "Both"))
            ),
            column(width = 4,
                actionButton(inputId = ns("featProcBut"), label = "Process Features:")
            )
        ),
        conditionalPanel("input.omic == 'Metabolomics'",
            fluidRow(
                column(width = 4,
                    switchInput(ns("ism0Filt"), label = "Filter non M0 features"),
                    switchInput(ns("annoFilt"), label = "Filter non annotated features")
                ),
                column(width = 4,
                    switchInput(inputId = ns("blankFilt"), label = "Blank Filter"),
                    switchInput(inputId = ns("cvFilt"), label = "Filter by variation coefficient"),
                    switchInput(inputId = ns("naFilt"), label = "Filter by NA ratio"),
                ),
                column(width = 4,
                    switchInput(inputId = ns("metvarFilt"), label = "Variance Filter"),
                    switchInput(inputId = ns("intFilt"), label = "Intensity Filter"),
                    switchInput(inputId = ns("sampleFilt"), label = "Sample Filter")
                )
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
                    conditionalPanel("input.blankFilt",
                        numericInput(inputId = ns("blankRatio"), label = "Select a fold-change between sample and blank:",
                                    value = 1, min = 0, max = 100, step = 0.1)

                    , ns = ns),
                    conditionalPanel("input.cvFilt",
                        numericInput(inputId = ns("maxRDS"),
                                    label = "Select a maximum QC RDS % threshold::",
                                    value = 20, min = 0, max = 100, step = 1)
                    , ns = ns),
                    conditionalPanel("input.naFilt",
                        numericInput(inputId = ns("nathr"), label = "NA ratio threhshold",
                                    value = 0.4, min = 0, max = 1, step = 0.01),
                        selectInput(inputId = ns("naMethod"), label = "Select a NA ratio method:",
                                    choices = c("WithinQC" = "QC", "WithinSampleClass" = "within",
                                                "AcrossAllSamples" = "across"))
                    , ns = ns)
                ),
                column(width = 4,
                    conditionalPanel("input.metvarFilt",
                        selectInput(inputId = ns("varfun"), label = "Select a variance function:",
                                    choices = c("IQR", "CV", "SD")),
                        switchInput(inputId = ns("quant"), label = "Filter by quantile"),
                        numericInput(inputId = ns("thr"), label = "Select variance threshold",
                                    value = 0.6, min = 0, max = 1000, step = 0.01),
                        selectInput(inputId = ns("groupName"),
                                    label = "Select a specific group to compute the variances (samples for example):",
                                    choices = "")
                    , ns = ns),
                    conditionalPanel("input.intFilt",
                        numericInput(inputId = ns("intthr"), label = "Select a intensity threhshold (absolute):",
                                    min = 0, value = 100, max = 10000000, step = 1)
                    , ns = ns),
                    conditionalPanel("input.sampleFilt",
                        numericInput(inputId = ns("maxNA"), label = "Select a missing value threshold:",
                                    min = 0, max = 1, value = 0.4, step = 0.01),
                        selectInput(inputId = ns("filtName"),
                                    label = "Select an optional group variable name to filter all it's samples:",
                                    choices = "")
                    , ns = ns)
                )
            )
        , ns = ns),
        conditionalPanel("input.omic == 'Transcriptomics'",
            fluidRow(
                column(width = 4,
                    switchInput(inputId = ns("entrez"), label = "Remove non entrez genes and duplicated")
                ),
                column(width = 4,
                    switchInput(inputId = ns("varfilt"), label = "Variance Filter")
                )
            ),
            fluidRow(
                conditionalPanel("input.varfilt",
                    column(width = 4,
                        selectInput(inputId = ns("varfun"), label = "Select the variance function to apply:",
                                    choices = c("IQR", "CV", "SD")),
                        numericInput(inputId = ns("quantvar"), label = "Select the quantile to eliminate:",
                                    value = 0.6, min = 0, max = 1, step = 0.01)
                    )
                ,ns = ns),
            )
        , ns = ns),
        conditionalPanel("input.omic == 'Both'",
            fluidRow(
                column(width = 4,
                    switchInput(inputId = ns("compFilter"), label = "Remove non-significant features")
                ),
                column(width = 4,
                    switchInput(inputId = ns("bioFilter"), label = "Filter with biosign model")
                )
            ),
            fluidRow(
                conditionalPanel("input.compFilter",
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
                , ns = ns),
                conditionalPanel("input.bioFilter",
                    column(width = 4,
                        selectInput(inputId = ns("bioSignMod"),
                                    label = "Select a biosign model to use for filtering:",
                                    choices = ""),
                        selectInput(inputId = ns("minScore"),
                                    label = "Select the minimum biosign score",
                                    choices = c("S","A", "B","E")),
                        selectInput(inputId = ns("bioSigMod"),
                                    label = "Select which model scores to use:",
                                    multiple = TRUE, choices = c("svm", "random forest", "pls-da"))
                    )
                ,ns = ns)
            )
        , ns = ns)
    )
}
