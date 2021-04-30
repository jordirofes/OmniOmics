processUI <- function(id){
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
            column(width = 4,
                actionButton(inputId = ns("proc"), label = "Process data:")
            )
        ),
        conditionalPanel("input.omic == 'Metabolomics'",
            fluidRow(
                column(width = 4,
                    switchInput(inputId = ns("rtcorrect"), label = "RTcorrection", value = TRUE),
                    switchInput(inputId = ns("peakmerge"), label = "Peak Merging", value = TRUE),
                    switchInput(inputId = ns("peakfill"), label = "Peak filling", value = TRUE)
                    ),
                column(width = 4,
                    radioButtons(inputId = ns("annotation"),
                                choiceNames = c("No annotation", "CAMERA", "cliqueMS"),
                                choiceValues = c("none", "camera", "cliqueMS"),
                                label = "Choose an annotation step:")
                )
            ),
            fluidRow(
                column(width = 3,
                    numericInput(inputId = ns("ppm"), label = "Allowed ppm error:",
                                value = 10, min = 1, max = 300, step = 1),
                    numericInput(inputId = ns("pwlb"), label = "Peak width lower bound:",
                                value = 10, min = 1, max = 200, step = 1),
                    numericInput(inputId = ns("pwup"), label = "Peak width upper bound:",
                                value = 120, min = 1, max = 200, step = 1),
                    numericInput(inputId = ns("noise"), label = "Noise param:",
                                value = 0, min = 0, max = 10000, step = 1),
                    numericInput(inputId = ns("snr"), label = "Minimum signal to noise ratio (SNR):",
                                value = 0, min = 0, max = 100000, step = 1)
                ),
                column(width = 3,
                    radioButtons(inputId = ns("phenoVar"), label = "Grouping Variable:",
                                choices = ""),
                    numericInput(inputId = ns("minfrac"), label = "Grouping minimun fraction of peaks:",
                                value = 0.4, min = 0, max = 1),
                    numericInput(inputId = ns("binwidth"), label = "Grouping binwidth:",
                                value = 1, min = 1, max = 100)
                ),
                column(width = 3,
                    conditionalPanel("input.peakmerge",
                        numericInput(inputId = ns("binsize"), label = "Peak merge bin size:",
                                    value = 1, min = 1, max = 100)
                    ,ns = ns),
                    conditionalPanel("input.rtcorrect",
                        numericInput(inputId = ns("expandrt"), label = "RTcorrection bindwith (Obiwarp):",
                                    value = 4, min = 0, max = 10, step = 1)
                    ,ns = ns),
                    conditionalPanel("input.annotation != 'none'",
                        radioButtons(inputId = ns("polarity"), label = "Polarity:",
                                    choices = c("positive", "negative"))
                    ,ns = ns),
                    conditionalPanel("input.annotation == 'cliqueMS'",
                        selectInput(inputId = ns("cliqueSample"),
                                    label = "Select which sample will be used for the cliqueMS annotation:",
                                    multiple = FALSE, choices = "")
                    ,ns = ns)
                )
            )
        , ns = ns),
        conditionalPanel("input.omic == 'Transcriptomics'",
            fluidRow(
                column(width = 4,
                    selectInput(inputId = "database", label = "Select an annotation gene database: (requires installed package)",
                                choices = "Default")
                )
            )
        , ns = ns)
    )
}
