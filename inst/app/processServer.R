processServer <- function(id, objectList){
    moduleServer(
        id,
        function(input, output, session){
            observeEvent({
                objectList$objects},{
                    objectNames <- 1:length(objectList$objects)
                    names(objectNames) <- names(objectList$objects)
                    updateSelectInput(inputId = "object",  choices = objectNames)
            }, ignoreInit = TRUE, ignoreNULL = TRUE, )

            observeEvent({input$object},{
                pheno_names <- colnames(extractPhenoData(objectList$objects[[as.numeric(input$object)]]))
                updateRadioButtons(inputId = "phenoVar", choices = pheno_names)
                if(class(objectList$objects[[input$object]]) == "onDiskMSnExp"){
                    updateSelectInput(inputId = "cliqueSample",
                                    choices = basename(fileNames(objectList$objects[[input$object]])))
                }
            }, ignoreInit = TRUE, ignoreNULL = TRUE)

            returnData <- reactiveValues(object = NULL, objectNames = NULL, trigger = 0, elimParentObj = FALSE)

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
