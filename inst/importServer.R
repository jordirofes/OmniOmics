importServer <- function(id){
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
                            )
                        ),
                        fluidRow(
                            conditionalPanel("input.omic == 'Metabolomics'",
                                             column(width = 4,
                                                    radioButtons(inputId = ns("phenoVar"), label = "Order Variable:",
                                                                 choices = "")
                                            ), ns = ns
                            ),
                            conditionalPanel("input.omic == 'Transcriptomics'",
                                             column(width = 4,
                                                    textInput(inputId = ns("geoData"), label = "GEOdataset valid entry:",
                                                            placeholder = "Ex. GSE46687")
                                            ), ns = ns
                            ),
                            column(width = 8,
                                    verbatimTextOutput(outputId = ns("path_names"), placeholder = TRUE)
                                )
                        )
                    )
                }, )
                updateRadioButtons(inputId = "phenoVar", choices = colnames(phenoData()))
            }, ignoreNULL = TRUE, ignoreInit = TRUE)
            loadedData <- reactive({
                validate(need(input$dt, message = FALSE))
                input$dt
            })

            shinyFileChoose(input, "dt_path", root=getVolumes(), filetypes=c("", "mzXML", "mzML"))

            data_paths <- reactive({
                if(is.integer(input$dt_path)){
                    "No files selected"
                } else{
                    parseFilePaths(roots = getVolumes(), selection = input$dt_path)$datapath
                }
            })
            output$path_names <- renderPrint({data_paths()})

            dtObject <- reactive({
                input$load_data
                isol_data_path <- isolate(data_paths())
                ext <- tools::file_ext(isol_data_path)
                ext <- tolower(ext)
                validate(need(isol_data_path, message = FALSE))
                if(isolate(input$omic) == "Metabolomics"){
                    validate(need(ext == "mzml" | ext == "mzxml", message = "Input a valid metabolomics data file"))
                    metaboImport(filedir = isol_data_path,
                                phenodata = isolate(phenoData()),
                                injectionvar = isolate({input$phenoVar}))
                } else if(isolate(input$omic) == "Transcriptomics"){
                    validate(need(ext == "cel", message = "", message = "Input a valid transcriptomic data file"))
                    transcriImport(datapath = isol_data_path,
                                phenodata = phenoData(), geoid = input$geoData)
                }
            })
            browser()
            observeEvent(eventExpr = dtObject(), {
                if(!is.null(dtObject)){
                    sendSweetAlert(title = "Data Loading",
                                text = "Your data was loaded successfully!",
                                type = "success", session = session)
                }
            }, ignoreInit = TRUE)
        }
    )
}
