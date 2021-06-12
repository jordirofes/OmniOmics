machineLearnServer <- function(id, objectList){
    moduleServer(id,
        function(input, output, session){
            returnData <- reactiveValues(object = NULL, objectNames = NULL, trigger = 0, elimParentObj = FALSE, elimNum = NULL)
            observeEvent({
                    objectList$objects},{
                        objectNames <- 1:length(objectList$objects)
                        names(objectNames) <- names(objectList$objects)
                        classDt <- lapply(objectList$objects, class)
                        objectNames2 <- objectNames[which(classDt == "SummarizedExperiment" | classDt == "ExpressionSet")]
                        objectNames3 <- objectNames[which(classDt == "train" | classDt == "biosign")]

                        updateSelectInput(inputId = "object",  choices = objectNames2)
                        updateSelectInput(inputId = "trainSet",  choices = objectNames2)
                        updateSelectInput(inputId = "trainSet2",  choices = objectNames2)
                        updateSelectInput(inputId = "finMod",  choices = objectNames3)
                        updateSelectInput(inputId = "testSet",  choices = objectNames2)
                }, ignoreInit = TRUE, ignoreNULL = TRUE)
            observeEvent({input$object},{
                obj <- objectList$objects[[as.numeric(input$object)]]
                pheno_names <- colnames(extractPhenoData(obj))
                updateSelectInput(inputId = "groupVar", choices = pheno_names)
            }, ignoreInit = TRUE, ignoreNULL = TRUE)

            observeEvent({input$testSet},{
                obj <- objectList$objects[[as.numeric(input$testSet)]]
                pheno_names <- colnames(extractPhenoData(obj))
                updateSelectInput(inputId = "groupVar2", choices = pheno_names)
            }, ignoreInit = TRUE, ignoreNULL = TRUE)

            observeEvent(input$groupVar2,{
                obj <- objectList$objects[[as.numeric(input$testSet)]]
                varNames <- extractPhenoData(obj)[[input$groupVar2]]
                updateSelectInput(inputId = "posClass", choices = as.character(unique(varNames)))
            }, ignoreInit = TRUE, ignoreNULL = TRUE)

            observeEvent(input$partition,{

                obj <- objectList$objects[[as.numeric(input$object)]]

                set.seed(input$randSeed)
                train <- sample(1:ncol(obj), size = round(ncol(obj)*input$trainSize), replace = FALSE)
                returnData$object <- list(obj[,train], obj[,-train])
                returnData$objectNames <- paste0(names(objectList$objects)[as.numeric(input$object)], c("_train_", "_test_"), c(input$trainSize, 1 - input$trainSize))
                returnData$trigger <- returnData$trigger + 1
            }, ignoreInit = TRUE, ignoreNULL = TRUE)

            observeEvent(input$biosignMod,{
                obj <- objectList$objects[[as.numeric(input$trainSet)]]
                varNames <- extractPhenoData(obj)[[input$groupVar2]]
                # validate(need(length(unique(varNames)) == 2, message = "Variable must have two levels"))
                returnData$object <- featureSign(features = obj, groupvar = input$groupVar)
                returnData$objectNames <- paste0(names(objectList$objects)[as.numeric(input$trainSet)], "_sign")
                returnData$trigger <- returnData$triger + 1
            }, ignoreInit = TRUE, ignoreNULL = TRUE)

            observeEvent(input$trainMod,{

                obj <- objectList$objects[[as.numeric(input$trainSet)]]

                returnData$object <- mlFit(dt = obj, groupvar = input$groupVar,
                                        method = input$mod,metric = input$metric,
                                        cvmethod = input$cvmeth,kfolds = input$fold,
                                        ntimes = input$partitions, preproc = input$prepro,
                                        tlength = input$tuneLength)
                returnData$objectNames <- paste0(names(objectList$objects)[as.numeric(input$trainSet)], "_",input$mod)
                returnData$trigger <- returnData$trigger + 1

            }, ignoreInit = TRUE, ignoreNULL = TRUE)

            observeEvent(input$plot,{
                obj <- objectList$objects[[as.numeric(input$testSet)]]
                mod <- objectList$objects[[as.numeric(input$finMod)]]
                output$cm <- renderPrint({
                    mlPredictCM(mlmod = mod, newdt = obj, groupvar = input$groupVar2,
                                posclass = input$posClass)

                })
                output$roc <- renderPlot({
                    mlPredictROC(mlmod = mod, newdt = obj, groupvar = input$groupVar2,
                                posclass = input$posClass)
                })
            }, ignoreInit = TRUE, ignoreNULL = TRUE)
            return(returnData)
        }
    )
}
