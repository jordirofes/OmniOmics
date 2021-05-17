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
                updateSelectInput(inputId = "phenoVar", choices = pheno_names)
                if(class(obj) == "OnDiskMSnExp"){
                    dt_files <- 1:length(basename(fileNames(obj)))
                    names(dt_files) <- basename(fileNames(obj))
                    updateSelectInput(inputId = "cliqueSample",
                                    choices = dt_files)
                }
                if(class(obj) == "GeneFeatureSet"){
                    packages <- installed.packages()[,1]
                    annotPack <- packages[grep(".db$", packages)]

                    validate(need(length(annotPack) != 0, message = "No annotation packages detected"))

                    selectedAnnotPack <- annotPack[1]

                    if(annotation(obj) %in% annotPack){
                        selectedAnnotPack <- annotation(obj)
                    }
                    updateSelectInput(inputId = "database", choices = annotPack, selected = selectedAnnotPack)
                }
            }, ignoreInit = TRUE, ignoreNULL = TRUE)

            observeEvent({input$objectDt}, {
                obj <- objectList$objects[[as.numeric(input$objectDt)]]
                validate(need(class(obj) == "OnDiskMSnExp" | class(obj) == "MSnExp", message = ""))
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
                obj <- objectList$objects[[as.numeric(input$objectDt)]]
                if(class(obj) == "OnDiskMSnExp" | class(obj) == "MSnExp" | class(obj) == "XCMSnExp"){
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
                                                cliqsamp = as.numeric(input$cliqueSample), mergepeaks = input$peakmerge,
                                                rtadjust = input$rtcorrect, group = TRUE, fill = input$peakfill)
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
