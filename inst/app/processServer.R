processServer <- function(id, objectList){
    moduleServer(
        id,
        function(input, output, session){
            observeEvent({
                objectList$objects},{
                    objectNames <- 1:length(objectList$objects)
                    names(objectNames) <- names(objectList$objects)
                    classDt <- lapply(objectList$objects, class)
                    objectNames1 <- objectNames[which(classDt == "OnDiskMSnExp" | classDt == "GeneFeatureSet" | classDt == "MSnExp")]
                    objectNames2 <- objectNames[which(classDt == "OnDiskMSnExp" | classDt == "GeneFeatureSet" | classDt == "MSnExp" | classDt == "ExpressionSet" | classDt == "XCMSnExp")]
                    updateSelectInput(inputId = "object",  choices = objectNames1)
                    updateSelectInput(inputId = "objectDt",  choices = objectNames2)
            }, ignoreInit = TRUE, ignoreNULL = TRUE, )

            observeEvent({input$object},{
                obj <- objectList$objects[[as.numeric(input$object)]]
                pheno_names <- colnames(extractPhenoData(obj))
                updateRadioButtons(inputId = "phenoVar", choices = pheno_names)
                if(class(objectList$objects[[input$object]]) == "onDiskMSnExp"){
                    updateSelectInput(inputId = "cliqueSample",
                                    choices = basename(fileNames(obj)))
                }
            }, ignoreInit = TRUE, ignoreNULL = TRUE)

            observeEvent({input$objectDt}, {
                obj <- objectList$objects[[as.numeric(input$objectDt)]]
                pheno_names <- colnames(extractPhenoData(obj))
                file_names <- 1:length(fileNames(obj))
                names(file_names) <- basename(fileNames(obj))
                updateSelectInput(inputId = "groupVar", choices = pheno_names)
                # updateSelectizeInput(inputId = "files", choices = file_names)
            }, ignoreInit = TRUE, ignoreNULL = TRUE)

            observeEvent({input$groupVar}, {
                obj <- objectList$objects[[as.numeric(input$objectDt)]]
                varNames <- extractPhenoData(obj)[[input$groupVar]]
                updateSelectInput(inputId = "groupFilt", choices = c("none", varNames))
            }, ignoreInit = TRUE, ignoreNULL = TRUE)

            observeEvent({input$objectDt}, {
                browser()
                obj <- objectList$objects[[as.numeric(input$objectDt)]]
                if(class(obj) == "OnDiskMSnEx" | class(obj) == "MSnEx" | class(obj) == "XCMSnExp"){
                    output$ggtic <- renderPlotly({
                        ggTicQuality(object = obj,
                                pheno_var = input$groupVar, pheno_filter = input$groupFilt,
                                violin = input$violin, order = input$order)
                    })
                } else if(class(obj) == "GeneFeatureSet" | class(obj) == "ExpressionSet"){
                    output$ggdistr <- renderPlotly({
                        ggExprDistrPlot(object = obj, groupvar = input$groupVar,
                                        violin = input$violin)
                    })
                    output$ggdens <- renderPlotly({
                        ggDensityPlot(object = obj, groupvar = input$groupVar)
                    })
                }
            }, ignoreInit = TRUE, ignoreNULL = TRUE)

            returnData <- reactiveValues(object = NULL, objectNames = NULL, trigger = 0)

            observeEvent({input$proc}, {
                validate(need(input$object, message = FALSE))
                objectName <- names(objectList$objects)[as.numeric(input$object)]
                if(input$omic == "Metabolomics"){
                    returnData$object <- metaboProc(object = objectList$objects[[as.numeric(input$object)]],
                                                polarity = input$polarity, groupvar = input$phenoVar,
                                                peakwidth = c(input$pwlb, input$pwup), noise = input$noise,
                                                snthresh = input$snr, ppm = input$ppm, expandrt = input$expandrt,
                                                binsize = input$binsize, minFraction = input$minfrac,
                                                bw = input$binwidth, annotation = input$annotation,
                                                cliqsamp = input$cliqueSample, mergepeaks = input$peakmerge,
                                                rtadjust = input$rtcorrect, group = TRUE, fill = input$peakfill,
                                                summ = TRUE)
                    returnData$objectNames <- paste0(objectName, "_processed")
                    if(input$annotation == "cliqueMS"){
                        returnData$objectNames <- c(returnData$objectNames, paste0(objectName, "_cliqueAnObj"))
                    }
                    if(input$annotation == "camera"){
                        returnData$objectNames <- c(returnData$objectNames, paste0(objectName, "_cameraAnObj"))
                    }
                    returnData$objectNames <- c(returnData$objectNames, paste0(objectName, "_summarizedExp"))
                }
                if(input$omic == "Transcriptomics"){
                    returnData$object <- procTranscript(object = objectList$objects[[as.numeric(input$object)]],
                                                        annotationTable = input$database)
                    returnData$objectNames <- paste0(objectName, "_processed")
                }
                sendSweetAlert(title = "Data Loading",
                                text = "Your data was processed successfully!",
                                type = "success", session = session)
                returnData$trigger <- returnData$trigger + 1
            })
            return(returnData)
        }
    )
}
