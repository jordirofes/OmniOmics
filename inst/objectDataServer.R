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
            }, ignoreInit = TRUE, ignoreNULL = TRUE, )
        }
    )
}
