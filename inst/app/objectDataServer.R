objectDataServer <- function(id, objectList){
    moduleServer(
        id,
        function(input, output, session){
            output$dt_info <- renderPrint({
                objectList$objects[[as.numeric(input$object)]]
            })
            observeEvent({
                objectList$objects},{
                    objectNames <- 1:length(objectList$objects)
                    names(objectNames) <- names(objectList$objects)
                    updateSelectInput(inputId = "object",  choices = objectNames)
            }, ignoreInit = TRUE, ignoreNULL = TRUE)

            observeEvent(input$object,{
                obj <- objectList$objects[[as.numeric(input$object)]]
                if(class(obj) == "SummarizedExperiment"){
                    dt_info <- as.data.frame(rowData(obj))
                    col_num <- which(colnames(dt_info) %in%  c("peakidx","annotation", "isotope"))
                    feat_dt <- as.matrix(dt_info[, -col_num])
                    output$extra_dt <- DT::renderDataTable(expr = {
                        feat_dt
                    }, options = list(pageLength = 10, scrollX = TRUE))
                }
                if(class(obj) == "XCMSnExp"){
                    output$extra_dt <- DT::renderDataTable(expr = {
                        chromPeaks(obj)
                    }, options = list(pageLength = 10, scrollX = TRUE))

                }
            })

        }
    )
}
