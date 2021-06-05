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
                    objectNames3 <- objectNames[which(classDt == "SummarizedExperiment" | classDt == "ExpressionSet")]
                    objectNames2 <- objectNames[which(classDt == "SummarizedExperiment" | classDt == "ExpressionSet" | classDt == "data.frame")]
                    objectNames1 <- objectNames[which(classDt == "data.frame")]
                    updateSelectInput(inputId = "compTable",  choices = objectNames1)
                    updateSelectInput(inputId = "object",  choices = objectNames3)
                    updateSelectInput(inputId = "object2",  choices = objectNames3)
                    updateSelectInput(inputId = "object3", choices = objectNames2)
            }, ignoreInit = TRUE, ignoreNULL = TRUE)

            observeEvent({input$object},{
                obj <- objectList$objects[[as.numeric(input$object)]]
                pheno_names <- colnames(extractPhenoData(obj))
                feat_names <- 1:nrow(extractData(obj))
                names(feat_names) <- rownames(extractData(obj))
                updateSelectizeInput(inputId = "feature", choices = feat_names)
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
                }, options = list(pageLength = 10, scrollX = TRUE))
            }, ignoreInit = TRUE, ignoreNULL = TRUE)

            observeEvent(input$contMatCalc,{
                obj <- objectList$objects[[as.numeric(input$object)]]
                modMat <- model.mat(obj, phenovar = input$groupVar2)

                contMat <- do.call(makeContrasts, list(input$cont1, input$cont2, input$cont3, levels = modMat))
                output$contMat <- renderDataTable({
                    contMat
                }, options = list(pageLength = 10, scrollX = TRUE))
            }, ignoreInit = TRUE, ignoreNULL = TRUE)

            observeEvent(input$compTable, {
                validate(need(input$compTable, message = ""))
                obj <- objectList$objects[[as.numeric(input$compTable)]]
                output$volcano <- renderPlotly({
                    groupFeatureVolcano(comptable = obj, adj.pvalue = input$adjpval)
                })
            }, ignoreNULL = TRUE, ignoreInit = TRUE)
            observeEvent(input$object2, {
                validate(need(input$object2, message = ""))
                obj <- objectList$objects[[as.numeric(input$object2)]]
                output$featComp <- renderPlot({
                    compPlot(obj, subset = as.numeric(input$feature), groupvar = input$groupVar3)
                })
            }, ignoreNULL = TRUE, ignoreInit = TRUE)
            observeEvent(input$compare,{
                obj <- objectList$objects[[as.numeric(input$object)]]
                if(input$omic == "Metabolomics"){
                    returnData$object <- groupComp(features = obj, groupvar = input$groupVar,
                                                elimlab = input$elimVar, test = input$compTest,
                                                adj.method = input$adjMethod, paired = input$paired,
                                                var.equal = input$eqVar)

                } else if(input$omic == "Transcriptomics"){
                    modMat <- model.mat(obj, phenovar = input$groupVar2)
                    contMat <- do.call(makeContrasts, list(input$cont1, input$cont2, input$cont3, levels = modMat))
                    returnData$object <- groupFeatureComp(features = obj, modelMatrix = modMat,
                                                    contrastMat = contMat, adjpval = input$adjMethod2)
                }
                gname <- lapply(returnData$object, function(x){ return(unique(x$compname))})
                returnData$objectNames <- paste0(names(objectList$objects)[as.numeric(input$object)],"_CT_", gname)
                returnData$trigger <- returnData$trigger + 1
            }, ignoreInit = TRUE, ignoreNULL = TRUE)

            observeEvent({input$object3
                input$omic2}, {
                validate(need(input$object3, message = ""))
                obj <- objectList$objects[[as.numeric(input$object3)]]
                validate(need(obj, message = ""))
                packages <- installed.packages()[,1]
                annotPack <- packages[grep(".db$", packages)]

                validate(need(length(annotPack) != 0, message = "No annotation packages detected"))
                selectedAnnotPack <- annotPack[1]
                if(class(obj) == "ExpressionSet"){
                    if(annotation(obj) %in% annotPack){
                        selectedAnnotPack <- annotation(obj)
                    }
                }
                if(input$omic2 == "Transcriptomics"){
                    updateSelectInput(inputId = "annotPack", choices = annotPack, selected = selectedAnnotPack)
                }
            }, ignoreNULL = TRUE, ignoreInit = TRUE)

            metFile <- reactive({
                validate(need(input$pd, message = FALSE))
                input$metList
            })
            metabData <- reactive({
                ext <- tools::file_ext(metFile()$datapath)
                req(metFile())
                validate(need(ext == "csv" | ext == "xls" | ext == "xlsx" | ext == "txt" | ext == "tsv",
                            message = "Please upload a valid metabolite data file"))
                validate(need(input$metList, message = FALSE))
                read.csv(metFile()$datapath)
            })

            observeEvent(input$annotate, {
                validate(need(input$object3, message = ""))
                validate(need(input$annotPack, message = ""))
                obj <- objectList$objects[[as.numeric(input$object3)]]
                do.call(library, list(input$annotPack))
                returnData$object <- annotateData(features = obj, anotpackage = input$annotPack, metabList = metabData())
                returnData$objectNames <- paste0(names(objectList$objects)[as.numeric(input$object3)], "_annotated")
                returnData$trigger <- returnData$trigger + 1

            }, ignoreNULL = TRUE, ignoreInit = TRUE)
            return(returnData)
        }
    )
}
