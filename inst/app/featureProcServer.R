featureProcServer <- function(id, objectList){
    moduleServer(
        id,
        function(input, output, session){
            returnData <- reactiveValues(object = NULL, objectNames = NULL, trigger = 0, elimParentObj = FALSE, elimNum = NULL)
            observeEvent({
                objectList$objects},{
                    objectNames <- 1:length(objectList$objects)
                    names(objectNames) <- names(objectList$objects)
                    updateSelectInput(inputId = "object",  choices = objectNames)
            }, ignoreInit = TRUE, ignoreNULL = TRUE)
            observeEvent({input$object},{
                obj <- objectList$objects[[as.numeric(input$object)]]
                pheno_names <- colnames(extractPhenoData(obj))
                updateRadioButtons(inputId = "groupVar", choices = pheno_names)
                updateRadioButtons(inputId = "groupVar2", choices = pheno_names)
            }, ignoreInit = TRUE, ignoreNULL = TRUE)
            observeEvent(input$groupVar,{
                obj <- objectList$objects[[as.numeric(input$object)]]
                varNames <- extractPhenoData(obj)[[input$groupVar]]
                updateSelectInput(inputId = "qcname", choices = as.character(unique(varNames)))
            }, ignoreInit = TRUE, ignoreNULL = TRUE)
        }
    )
}
