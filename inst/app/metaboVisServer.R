metaboVisServer <- function(id, objectList){
    moduleServer(id,
        function(input, output, session){
            observeEvent({
                objectList$objects},{
                    objectNames <- 1:length(objectList$objects)
                    names(objectNames) <- names(objectList$objects)
                    classDt <- lapply(objectList$objects, class)
                    objectNames <- objectNames[which(classDt == "OnDiskMSnExp" | classDt == "MSnExp")]
                    updateSelectInput(inputId = "object",  choices = objectNames)
            }, ignoreInit = TRUE, ignoreNULL = TRUE)

            observeEvent({input$object},{
                obj <- objectList$objects[[as.numeric(input$object)]]
                pheno_names <- colnames(extractPhenoData(obj))
                updateSelectInput(inputId = "groupVar", choices = pheno_names)

                file_names <- 1:length(basename(fileNames(obj)))
                names(file_names) <- basename(fileNames(obj))
                updateSelectInput(inputId = "files", choices = file_names)
                rt_max <- round(max(obj@featureData@data$retentionTime), 5)
                rt_min <- round(min(obj@featureData@data$retentionTime), 5)
                updateSliderInput(inputId = "rt", min = rt_min, max = rt_max, value = c(rt_min, rt_max))
                mz_min <- round(min(obj@featureData@data$lowMZ), 5)
                mz_max <- round(max(obj@featureData@data$highMZ), 5)
                updateNumericInput(inputId = "mz2", min = mz_min, max = mz_max, value = ((mz_min + mz_max)/2), step = 0.0001)
                updateSliderInput(inputId = "mz1", min = mz_min, max = mz_max, value = c(mz_min, mz_max), step = 0.0001)
            }, ignoreInit = TRUE, ignoreNULL = TRUE)

            observe({
                validate(need(input$object, message = ""))
                validate(need(input$files, message = ""))
                obj <- objectList$objects[[as.numeric(input$object)]]
                if(input$byrange){
                    mz <- input$mz1
                } else{
                    mz <- input$mz2
                }
                output$chromplot <- renderPlotly({
                    ggChromPlot(object = obj, filenum = as.numeric(input$files),
                                mz = mz, ppm = input$ppm, rtint = input$rt,
                                pheno_var = input$groupVar, chromtype = input$chromtype,
                                logscale = input$log)
                })
            })
        }
    )
}
