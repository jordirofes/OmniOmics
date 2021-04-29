processUI <- function(id){
    ns <- NS(id)
    tagList(
        fileInput(inputId = ns("pd"), label = "Load your Pheno Data:",
                buttonLabel = "Load", accept = c(".csv",".txt", ".xls",
                                                ".xlsx"))
    )
}
