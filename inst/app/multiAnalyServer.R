multiAnalyServer <- function(id, objectList){
    moduleServer(id,
        function(input, output, session){
            returnData <- reactiveValues(object = NULL, objectNames = NULL, trigger = 0, elimParentObj = FALSE, elimNum = NULL)
            observeEvent({
                    objectList$objects},{
                        objectNames <- 1:length(objectList$objects)
                        names(objectNames) <- names(objectList$objects)
                        classDt <- lapply(objectList$objects, class)
                        objectNames2 <- objectNames[which(classDt == "SummarizedExperiment" | classDt == "ExpressionSet" | classDt == "GeneFeatureSet")]
                        objectNames3 <- objectNames[which(classDt == "prcomp")]
                        updateSelectInput(inputId = "object",  choices = objectNames2)
                        updateSelectInput(inputId = "object2",  choices = objectNames2)
                        updateSelectInput(inputId = "pcaObj",  choices = objectNames3)
                }, ignoreInit = TRUE, ignoreNULL = TRUE)

            observeEvent({input$object2},{
                obj <- objectList$objects[[as.numeric(input$object2)]]
                pheno_names <- colnames(extractPhenoData(obj))
                updateSelectInput(inputId = "groupVar2", choices = pheno_names)
            }, ignoreInit = TRUE, ignoreNULL = TRUE)

            observeEvent({input$object},{
                obj <- objectList$objects[[as.numeric(input$object)]]
                pheno_names <- colnames(extractPhenoData(obj))
                updateSelectInput(inputId = "groupVar", choices = pheno_names)
            }, ignoreInit = TRUE, ignoreNULL = TRUE)

            observeEvent({input$pcaCalc}, {
                obj <- objectList$objects[[as.numeric(input$object)]]
                returnData$object <- multipca(obj, scale = input$scaling)
                sendSweetAlert(title = "Data processing",
                                text = "Your pca was calculated successfully!",
                                type = "success", session = session)
                returnData$objectNames <- paste0(names(objectList$objects)[as.numeric(input$object)], "_pca")
                returnData$trigger <- returnData$trigger + 1
            }, ignoreInit = TRUE, ignoreNULL = TRUE)

            observeEvent({
                input$pcaObj
                input$object
            },{
                obj <- objectList$objects[[as.numeric(input$object)]]
                pcaObj <- objectList$objects[[as.numeric(input$pcaObj)]]

                if(class(pcaObj) == "prcomp"){
                    updateSelectInput(inputId = "pc1", choices = 1:ncol(pcaObj$x))
                    updateSelectInput(inputId = "pc2", choices = 1:ncol(pcaObj$x), selected = 2)
                    updateSelectInput(inputId = "pcload", choices = 1:ncol(pcaObj$x))
                    output$pcaPlot <- renderPlotly(pcaPlot(pcdt = pcaObj, object = obj,
                                                    pc = c(as.numeric(input$pc1), as.numeric(input$pc2)),
                                                    groupvar = input$groupVar, interactive = TRUE,
                                                    bygroup = input$grouped))
                }
            }, ignoreInit = TRUE, ignoreNULL = TRUE)
            observeEvent({input$pcload}, {
                pcaObj <- objectList$objects[[as.numeric(input$pcaObj)]]
                output$loadPlot <- renderPlotly(loadingsPlot(pcdt = pcaObj, pc = as.numeric(input$pcload)))
            })

            observeEvent({
                input$object2
                input$groupVar2
            },{
                obj <- objectList$objects[[as.numeric(input$object2)]]
                if(class(obj) == "SummarizedExperiment" | class(obj) == "ExpressionSet"){
                    output$heatPlot <- renderPlot(heatmapPlot(obj, groupvar = input$groupVar2))
                }
            }, ignoreInit = TRUE, ignoreNULL = TRUE)

            return(returnData)
        }
    )
}
