groupCompServer <- function(id, objectList){
    moduleServer(
        id,
        function(input, output, session){
            returnData <- reactiveValues(object = NULL, objectNames = NULL, trigger = 0, elimParentObj = FALSE, elimNum = NULL)
            observeEvent({
                objectList$objects},{
                    objectNames <- 1:length(objectList$objects)
                    names(objectNames) <- names(objectList$objects)
                    classDt <- lapply(objectList$objects, class)
                    objectNames <- objectNames[which(classDt == "SummarizedExperiment" | classDt == "ExpressionSet")]
                    objectNames1 <- objectNames[which(classDt == "data.frame")]
                    updateSelectInput(inputId = "compTable",  choices = objectNames1)
                    updateSelectInput(inputId = "object",  choices = objectNames)
                    updateSelectInput(inputId = "object2",  choices = objectNames)
            }, ignoreInit = TRUE, ignoreNULL = TRUE)

            observeEvent({input$object},{
                obj <- objectList$objects[[as.numeric(input$object)]]
                pheno_names <- colnames(extractPhenoData(obj))
                feat_names <- 1:length(nrow(obj))
                names(feat_names) <- nrow(obj)
                updateSelectInput(inputId = "feature", choices = feat_names)
                updateSelectInput(inputId = "groupVar", choices = pheno_names)
                updateSelectInput(inputId = "groupVar2", choices = pheno_names)
                updateSelectInput(inputId = "groupVar3", choices = pheno_names)
            }, ignoreInit = TRUE, ignoreNULL = TRUE)

            observeEvent(input$groupVar,{
                obj <- objectList$objects[[as.numeric(input$object)]]
                varNames <- extractPhenoData(obj)[[input$groupVar]]
                updateSelectInput(inputId = "elimVar", choices = c("none", as.character(unique(varNames))))
            }, ignoreInit = TRUE, ignoreNULL = TRUE)

            observeEvent(input$groupVar2,{
                obj <- objectList$objects[[as.numeric(input$object)]]
                modMat <- model.mat(obj, phenovar = input$groupVar2)
                output$compMat <- renderDataTable({
                    modMat
                })
            }, ignoreInit = TRUE, ignoreNULL = TRUE)

            observeEvent(input$contMatCalc,{
                contMat <- makeContrasts(cont1, cont2, cont3, levels = modMat)
                output$contMat <- renderDataTable({
                    contMat
                }, options = list(pageLength = 10, scrollX = TRUE))
            }, ignoreInit = TRUE, ignoreNULL = TRUE)

            observeEvent(input$compTable, {
                obj <- objectList$objects[[as.numeric(input$compTable)]]
                output$volcano <- renderPlotly({
                    groupFeatureVolcano(comptable = obj, adj.pvalue = input$adjpval)
                })
            })
            observeEvent(input$object2, {
                obj <- objectList$objects[[as.numerci(input$object2)]]
                output$featComp <- renderPlot({
                    compPlot(obj, features = input$feature, groupvar = input$groupVar3)
                })
            })
            observeEvent(input$compare,{
                obj <- objectList$objects[[as.numeric(input$object)]]
                if(input$omic == "Metabolomics"){
                    returnData$object <- groupComp(features = obj, groupvar = input$groupVar,
                                                elimlab = input$elimVar, test = input$compTest,
                                                adj.method = input$adjMethod, paired = input$paired,
                                                var.equal = input$eqVar)


                } else if(input$omic == "Transcriptomics"){
                    returnData$object <- groupFeatureComp(features = obj, modelMatrix = modMat,
                                                    contrastMat = contMat, adjpval = input$adjMethod2)
                }
                gname <- lapply(returnData$object, function(x){ return(unique(x$compname))})
                returnData$objectNames <- paste0(names(objectList$objects)[as.numeric(input$object)],"_CT_", gname)
                returnData$trigger <- returnData$trigger + 1
            }, ignoreInit = TRUE, ignoreNULL = TRUE)
        }
    )
}
