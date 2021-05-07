multiAnalyUI <- function(id){
    ns <- NS(id)
    fluidPage(
        tabsetPanel(
            tabPanel(title = "PCA",
                fluidRow(
                    column(width = 4,
                        selectInput(inputId = ns("object"), label = "Select an object to process:",
                                    choices = ""),
                        selectInput(inputId = ns("pcaObj"), label = "Input the PCA object:",
                                    choices = "")
                    ),
                    column(width = 4,
                        switchInput(inputId = ns("scaling"), label = "PCA scaling")
                    ),
                    column(width = 4,
                        actionButton(inputId = ns("pcaCalc"), label = "Calculate PCA")
                    )
                ),
                fluidRow(
                    column(width = 4,
                        selectInput(inputId = ns("pc1"), label = "PCs to plot:", choices = ""),
                        selectInput(inputId = ns("pc2"), label = "", choices = ""),
                        switchInput(inputId = ns("grouped"), label = "By group")

                    ),
                    column(width = 4,
                        selectInput(inputId = ns("pcload"), label = "Loading PC:", choices = "")
                    )
                ),
                fluidRow(
                    column(width = 12,
                        selectInput(inputId = ns("groupVar"), label = "", choices = "")
                    )
                ),
                fluidRow(
                    plotlyOutput(outputId = ns("pcaPlot"))
                ),
                fluidRow(
                    plotlyOutput(outputId = ns("loadPlot"))
                )
            ),
            tabPanel(title = "Heatmap",
                fluidRow(
                    column(width = 4,
                        selectInput(inputId = ns("object2"), label = "Select an object to process:",
                                    choices = "")
                    ),
                    column(width = 4,
                        selectInput(inputId = ns("groupVar2"), label = "", choices = "")
                    )
                ),
                fluidRow(
                    column(width = 12,
                        plotOutput(outputId = ns("heatPlot"))
                    )
                ,)
            )
        ,type = "tabs")
    )
}
