metaboVistServer <- function(id, objectList){
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
                rt_max <- max(rtime(obj[[1]]))
                rt_min <- min(rtime(obj[[1]]))
                updateNumericInput(inputId = "rt", min = rt_min, max = rt_max, value = c(rt_min, rt_max))

            }, ignoreInit = TRUE, ignoreNULL = TRUE)

            observeEvent(input$object, {
                obj <- objectList$objects[[as.numeric(input$object)]]
                output$chromplot <- renderPlotly({
                    ggChromPlot(object = obj, filenum = input$files,
                                mz = input$mz, ppm = input$ppm, rtint = input$rt,
                                pheno_var = input$groupVar, chromtype = input$chromtype)
                })
            })
        }
    )
}
