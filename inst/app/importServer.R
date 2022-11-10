importServer <- function(id, objectList){
    moduleServer(
        id,
        function(input, output, session){
            phenoFile <- reactive({
                validate(need(input$pd, message = FALSE))
                input$pd
            })
            phenoData <- reactive({
                ext <- tools::file_ext(phenoFile()$datapath)
                req(phenoFile())
                validate(need(ext == "csv" | ext == "xls" | ext == "xlsx" | ext == "txt" | ext == "tsv",
                            message = "Please upload a valid pheno data file"))
                validate(need(input$pd, message = FALSE))
                phenoImport(phenoFile()$datapath, isolate({input$header}),
                        isolate({input$sep}))
            })

            output$phenoTable <- DT::renderDataTable(expr = {
                validate(need(input$pd, message = FALSE))
                phenoData()
            }, options = list(pageLength = 10, scrollX = TRUE))

            observeEvent(input$pd, {
                output$dt_upload_ui <- renderUI(expr = {
                    ns <- session$ns
                    tagList(
                        fluidRow(
                            column(width = 4,
                                   selectInput(inputId = ns("omic"), label = "Select your omic:",choices = c("Metabolomics", "Transcriptomics"))
                            ),
                            column(width = 4,
                                shinyFilesButton(ns("dt_path"),style = "margin-top: 25px;", label="File select", title="Please select a file", multiple=TRUE),
                                actionButton(ns("load_data"), label = "Load", style = "margin-top: 25px;")
                            ),
                            column(width = 4,
                                textInput(inputId = ns("fileName"), value = "experiment1", label = "Enter a name for the loaded object:"))
                        ),
                        fluidRow(
                            conditionalPanel("input.omic == 'Metabolomics'",
                                            column(width = 4,
                                                    radioButtons(inputId = ns("phenoVar"), label = "Order Variable:",
                                                                choices = "", selected = character(0))
                                            ), ns = ns
                            ),
                            conditionalPanel("input.omic == 'Transcriptomics'",
                                             column(width = 4,
                                                    textInput(inputId = ns("geoData"), label = "GEOdataset valid entry:",
                                                            placeholder = "Ex. GSE46687", value = NA)
                                            ), ns = ns
                            ),
                            column(width = 8,
                                    verbatimTextOutput(outputId = ns("path_names"), placeholder = TRUE)
                            )
                        )
                    )
                }, )
                updateRadioButtons(inputId = "phenoVar", choices = c("rownames", colnames(phenoData())))
            }, ignoreNULL = TRUE, ignoreInit = TRUE)

            loadedData <- reactive({
                validate(need(input$dt, message = FALSE))
                input$dt
            })

            shinyFileChoose(input, "dt_path", root=getVolumes(), filetypes=c("", "mzXML", "mzML", "CEL"))

            data_paths <- reactive({
                    parseFilePaths(roots = getVolumes(), selection = input$dt_path)$datapath
            })
            output$path_names <- renderPrint({
                if(is.integer(input$dt_path)){
                    "No files selected"
                } else{
                    data_paths()
                }
            })
            returnData <- reactiveValues(object = NULL, trigger = 0)
            observeEvent(input$load_data, {
                ext <- tools::file_ext(isolate(data_paths()))
                ext <- tolower(ext)
                if(isolate(input$omic) == "Metabolomics"){
                    validate(need(ext == "mzml" | ext == "mzxml", message = "Input a valid metabolomics data file"))
                    returnData$object <- metaboImport(filedir = isolate(data_paths()),
                                                phenodata = isolate(phenoData()),
                                                injectionvar = isolate({input$phenoVar}))
                } else if(isolate(input$omic) == "Transcriptomics"){
                    validate(need(ext == "cel", message = ""))
                    returnData$object <- transcriImport(datapath = isolate(data_paths()),
                                    phenodata = isolate(phenoData()),
                                    geoid = isolate(input$geoData))
                }
                sendSweetAlert(title = "Data Loading",
                                text = "Your data was loaded successfully!",
                                type = "success", session = session)
                returnData$trigger <- isolate(returnData$trigger) + 1
                returnData$objectName <- input$fileName
            }, ignoreNULL = TRUE, ignoreInit = TRUE, priority = 1)

            shinyFileChoose(input, "dt_path2", root=getVolumes(), filetypes=c("rds"))

            data_paths2 <- reactive({
                    parseFilePaths(roots = getVolumes(), selection = input$dt_path2)$datapath
            })
            output$path_names2 <- renderPrint({
                if(is.integer(input$dt_path2)){
                    "No files selected"
                } else{
                    data_paths2()
                }
            })
            observeEvent(input$load_data2, {
                ext <- tools::file_ext(isolate(data_paths2()))
                ext <- tolower(ext)
                if(ext == "rds"){
                    returnData$object <- readRDS(file = isolate(data_paths2()))
                }
                sendSweetAlert(title = "Data Loading",
                                text = "Your data was loaded successfully!",
                                type = "success", session = session)
                returnData$trigger <- isolate(returnData$trigger) + 1
                if(!is(returnData$object, class2 = "list")){
                    returnData$object <- list(returnData$object)
                }
                if(length(returnData$object) != 1){
                    returnData$objectName <- paste0(input$fileName2, 1:length(returnData$object))
                } else{
                    returnData$objectName <- input$fileName2
                }
            }, ignoreNULL = TRUE, ignoreInit = TRUE, priority = 1)

            observeEvent({
                objectList$objects},{
                    objectNames <- 1:length(objectList$objects)
                    names(objectNames) <- names(objectList$objects)
                    updateSelectInput(inputId = "object",  choices = objectNames)
            }, ignoreInit = TRUE, ignoreNULL = TRUE)

            shinyFileSave(input, "folder", roots= getVolumes(), filetypes=c(""))
            save_path <- reactive(parseSavePath(getVolumes(), input$folder)$datapath)

            output$dir_name <- renderPrint({
                if(length(save_path()) == 0){
                    "No files selected"
                } else{
                    save_path()
                }
            })

            observeEvent(input$save_data, {
                validate(need(input$folder, message = ""))
                validate(need(input$object, message = ""))
                obj <- objectList$objects[as.numeric(input$object)]
                saveRDS(object = obj, file = paste0(save_path(), ".rds"))
                sendSweetAlert(title = "Data Saving",
                                text = "Your data was saved successfully!",
                                type = "success", session = session)
            }, ignoreInit = TRUE, ignoreNULL = TRUE)
            return(returnData)
        }
    )
}
