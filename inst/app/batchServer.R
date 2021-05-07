batchServer <- function(id, objectList){
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
                    updateSelectInput(inputId = "object",  choices = objectNames)
                    updateSelectInput(inputId = "object2",  choices = objectNames)
            }, ignoreInit = TRUE, ignoreNULL = TRUE)

            observeEvent({input$object},{
                obj <- objectList$objects[[as.numeric(input$object)]]
                pheno_names <- colnames(extractPhenoData(obj))
                updateRadioButtons(inputId = "phenoVar", choices = pheno_names)
                updateRadioButtons(inputId = "ordVar", choices = c("Default", pheno_names))
                updateRadioButtons(inputId = "batchVar", choices = c("Default", pheno_names))
                updateRadioButtons(inputId = "phenoVar2", choices = pheno_names)
                updateRadioButtons(inputId = "groupVar", choices = pheno_names)
            }, ignoreInit = TRUE, ignoreNULL = TRUE)

            observeEvent(input$groupVar,{
                obj <- objectList$objects[[as.numeric(input$object)]]
                varNames <- extractPhenoData(obj)[[input$groupVar]]
                updateSelectInput(inputId = "qcname", choices = as.character(unique(varNames)))
            }, ignoreInit = TRUE, ignoreNULL = TRUE)

            observeEvent({input$object2},{
                obj <- objectList$objects[[as.numeric(input$object2)]]
                pheno_names <- colnames(extractPhenoData(obj))
                updateSelectInput(inputId = "groupVar2", choices = pheno_names)
                updateSelectInput(inputId = "phenoVar3", choices = pheno_names)
            }, ignoreInit = TRUE, ignoreNULL = TRUE)

            observeEvent(input$groupVar2,{
                obj <- objectList$objects[[as.numeric(input$object2)]]
                varNames <- extractPhenoData(obj)[[input$groupVar2]]
                updateSelectInput(inputId = "qcname2", choices = as.character(unique(varNames)))
            }, ignoreInit = TRUE, ignoreNULL = TRUE)

            observeEvent(input$groupVar2, {
                obj <- objectList$objects[[as.numeric(input$object2)]]
                if(class(obj) == "SummarizedExperiment"){
                    output$injplot <- renderPlotly({
                        featurebatchQc(features = obj, groupvar = input$groupVar2,
                                    qcname = input$qcname2, interactive = TRUE)
                    })
                } else if(class(obj) == "ExpressionSet"){
                    output$covplot <- renderPlotly({
                        featureBatchPVCA(features = obj, phenovars = input$phenoVar3,
                                        threshold = input$thr)
                    })
                }
            })

            observeEvent(input$batchBut,{
                obj <- objectList$objects[[as.numeric(input$object)]]
                if(input$batchNorm == "Covariate"){
                    meth <- "covnorm"
                }
                if(input$batchNorm == "Injection Order"){
                    meth <- "qcnorm"
                }
                if(input$ordVar == "Default"){
                    injOrder <- 1:nrow(obj)
                } else{
                    injOrder <- extractPhenoData(obj)[[input$ordVar]]
                }
                if(input$batchVar == "Default"){
                    batchOrd <- rep(1, nrow(obj))
                } else{
                    batchOrd <- extractPhenoData(obj)[[input$batchVar]]
                }
                returnData$object <- batchNormalization(features = obj, method = meth,
                                        injectionorder = injOrder, batchnum = batchOrd,
                                        groups = input$groupVar, qcname = input$qcname,
                                        covariate = input$phenoVar, covariate2 = input$phenoVar2)
                sendSweetAlert(title = "Data Loading",
                                text = "Your data was batch corrected successfully!",
                                type = "success", session = session)
                returnData$trigger <- returnData$trigger + 1
                returnData$objectNames <- paste0(objectNames, "_BC")
            }, ignoreInit = TRUE, ignoreNULL = TRUE)
            return(returnData)
        }
    )
}
