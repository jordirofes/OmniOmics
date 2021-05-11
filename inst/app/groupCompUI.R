groupCompUI <- function(id){
    ns <- NS(id)
    tabsetPanel(
        tabPanel(
            fluidPage(
                fluidRow(
                    column(width = 4,
                        selectInput(inputId = ns("object"), label = "Select an object to process:",
                                    choices = "")
                    ),
                    column(width = 4,
                        selectInput(inputId = ns("omic"), label = "Select group comparison methodology:",
                                    choices = c("Metabolomics", "Transcriptomics"))
                    ),
                    column(width = 4,
                        actionButton(inputId = ns("compare"), label = "Create comparison tables")
                    )
                ),
                conditionalPanel("input.omic == 'Metabolomics'",
                    fluidRow(
                        column(width = 4,
                            selectInput(inputId = ns("groupVar"), label = "Grouping variable to make comparisons:",
                                        choices = "")
                        ),
                        column(width = 4,
                            selectInput(inputId = ns("compTest"), label = "Select a test to compare groups:",
                                        choices = c("t-test", "mann-whitney")),
                            switchInput(inputId = ns("paired"), label = "Do a paired t-test"),
                            switchInput(inputId = ns("eqVar"), label = "Equal variance"),
                            selectInput(inputId = ns("adjMethod"), label = "Select a p-value adjust method:",
                                        choices = c("fdr", "bonferroni")),
                            selectInput(inputId = ns("elimVar"), label = "Eliminate one of the group labs: (optional)",
                                        choices = "")
                        )
                    )
                , ns = ns),
                conditionalPanel("input.omic == 'Transcriptomics'",
                    fluidRow(
                        column(width = 4,
                            selectInput(inputId = ns("groupVar2"), label = "GroupVar to create the comparison matrix:",
                                        choices = "")
                        ),
                        column(width = 4,
                            actionButton(ns("contMatCalc"), label = "Calculate contrast matrix"),
                            textInput(inputId = ns("cont1"), label = "Input a contrast between groups (use the names of the comparison matrix):",
                                    placeholder = c("GroupA - GroupB")),
                            textInput(inputId = ns("cont1"), label = "",
                                    placeholder = c("GroupC - GroupB")),
                            textInput(inputId = ns("cont1"), label = "",
                                    placeholder = c("(GroupA - GroupB) - (GroupC - Group)"))
                        ),
                        column(width = 4,
                            selectInput(inputId = ns("adjMethod2"), label = "Select a p-value adjust method:",
                                        choices = c("fdr", "bonferroni"))
                        )
                    ),
                    fluidRow(
                        column(width = 6,
                            DT::DTOutput(outputId = ns("compMat"))
                        ),
                        column(width = 6,
                            DT::DTOutput(outputId = ns("contMat"))
                        )
                    )
                ,ns = ns)
            )
        ),
        tabPanel(title = "Volcano Plot",
            fluidPage(
                fluidRow(
                    column(width = 4,
                        selectInput(inputId = ns("compTable"), label = "Select a comparison table:",
                                    choices = "")
                    ),
                    column(width = 4,
                        switchInput(inputId = ns("adjpval"), label = "Adjusted p-value")
                    )
                ),
                fluidRow(
                    column(width = 12,
                        plotlyOutput(outputId = ns("volcano"))
                    )
                )
            )
        ),
        tabPanel(title = "Features Comp",
            fluidPage(
                fluidRow(
                    column(width = 4,
                        selectInput(inputId = ns("object2"), label = "Select an object to process:",
                                    choices = "")
                    ),
                    column(width = 4,
                        selectInput(inputId = ns("groupVar3"), label = "Group variable:",
                                    choices = "")
                    ),
                    column(width = 4,
                        selectInput(inputId = ns("feature"), label = "Feature to compare:",
                                    choices = "")
                    )
                ),
                fluidRow(
                    column(width = 12,
                        plotOutput(outputId = ns("featComp"))
                    )
                )
            )
        )
    )

}
