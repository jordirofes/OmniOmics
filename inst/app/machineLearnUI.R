machineLearnUI <- function(id){
    ns <- NS(id)
    fluidPage(
        tabsetPanel(
            tabPanel(title = "Training",
                fluidRow(
                    column(width = 4,
                           selectInput(inputId = ns("object"), label = "Select an object to process:",
                                       choices = ""),
                           selectInput(inputId = ns("groupVar"), label = "Select the outcome variable",
                                       choices = "")
                    ),
                    column(width = 4,
                           numericInput(inputId = ns("trainSize"), label = "Train set size:",
                                        value = 0.7, min = 0, max = 1, step = 0.01),
                           numericInput(inputId = ns("randSeed"), label = "Train set size:",
                                        value = 123, min = 0, max = 1000000, step = 1)
                    ),
                    column(width = 4,
                           actionButton(inputId = ns("partition"), label = "Create partition")
                    )
                ),
                tabsetPanel(
                    tabPanel(title = "Biosign",
                        fluidRow(
                            column(width = 4,
                                selectInput(inputId = ns("trainSet"), label = "Select the training set:",
                                            choices = "")
                            ),
                            column(width = 4,
                                actionButton(inputId = ns("biosignMod"), label = "Create the biosign model")
                            )
                        ),
                        fluidRow(
                            column(width = 12,
                                plotOutput(outputId = ns("biosignPlot"))
                            )
                        )
                    ,),
                    tabPanel(title = "Caret",
                        fluidRow(
                            column(width = 4,
                                selectInput(inputId = ns("trainSet2"), label = "Select the training set:",
                                            choices = "")
                            ),
                            column(width = 4,
                                actionButton(inputId = ns("trainMod"), label = "Start training")
                            )
                        ),
                        fluidRow(
                            column(width = 4,
                                selectInput(inputId = ns("mod"), label = "Training model:",
                                            choices = names(getModelInfo())),
                                selectInput(inputId = ns("metric"), label = "Training metric:",
                                            choices = c("Accuracy", "Kappa")),
                                selectInput(inputId = ns("prepro"), label = "Pre-processing steps:",
                                            choices = c("center", "scale", "range", "knnImpute","BoxCox"), multiple = TRUE),
                                selectInput(inputId = ns("cvmeth"), label = "Cross-validation method:",
                                            choices = c("k-fold", "loo"))
                            ),
                            column(width = 4,
                                numericInput(inputId = ns("fold"), label = "K-folds:",
                                            value = 10, min = 1, max = 100, step = 1),
                                numericInput(inputId = ns("partitions"), label = "Number of partitions:",
                                            value = 5, min = 1, max = 100, step = 1),
                                numericInput(inputId = ns("tuneLength"), label = "Granularity of tuning parameters grid:",
                                            value = 10, min = 1, max = 100, step = 1)
                            )
                        ),
                    )
                )
            ),
            tabPanel(title = "Test & Plot",
                fluidRow(
                    column(width = 4,
                        selectInput(inputId = ns("finMod"), label = "Select the model to predict:",
                                    choices = ""),
                        selectInput(inputId = ns("groupVar2"), label = "Grouping Variable:",
                                    choices = "")
                    ),
                    column(width = 4,
                        selectInput(inputId = ns("testSet"), label = "Select the test set:",
                                    choices = ""),
                        selectInput(inputId = ns("posClass"), label = "Select the positive class:",
                                    choices = "")
                    ),
                    column(width = 4,
                        actionButton(inputId = ns("plot"), label = "Calculate results")
                    )
                ),
                fluidRow(
                    column(width = 6,
                        verbatimTextOutput(outputId = ns("cm"))
                    ),
                    column(width = 6,
                        plotOutput(outputId = ns("roc"))
                    )
                )
            )
        ),
    )
}
