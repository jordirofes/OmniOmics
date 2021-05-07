featureProcServer <- function(id, objectList){
    moduleServer(
        id,
        function(input, output, session){
            returnData <- reactiveValues(object = NULL, objectNames = NULL, trigger = 0, elimParentObj = FALSE, elimNum = NULL)
            observeEvent({
                objectList$objects},{
                    objectNames <- 1:length(objectList$objects)
                    names(objectNames) <- names(objectList$objects)
                    classDt <- lapply(objectList$objects, class)
                    objectNames <- objectNames[which(classDt == "SummarizedExperiment" | classDt == "ExpressionSet")]
                    updateSelectInput(inputId = "object",  choices = objectNames)
                    updateSelectInput(inputId = "object2",  choices = objectNames)
                    updateSelectInput(inputId = "bioSign",  choices = objectNames)
                    updateSelectInput(inputId = "compTable",  choices = objectNames)

            }, ignoreInit = TRUE, ignoreNULL = TRUE)
            observeEvent({input$object},{
                obj <- objectList$objects[[as.numeric(input$object)]]
                pheno_names <- colnames(extractPhenoData(obj))
                updateSelectInput(inputId = "groupVar", choices = pheno_names)
                updateSelectInput(inputId = "groupVar2", choices = pheno_names)
            }, ignoreInit = TRUE, ignoreNULL = TRUE)

            observeEvent(input$object2, {
                obj <- objectList$objects[[as.numeric(input$object2)]]
                if(input$varfunction == "IQR"){
                    varfun <- IQR
                } else if(input$varfunction == "CV"){
                    varfun <- function(x){sd(x)/mean(x)}
                } else if(input$varfunction == "SD"){
                    varfun <- sd
                }
                output$varplot <- renderPlotly({
                    featureVarPlot(features = obj, varfun = varfun)
                })
            })

            observeEvent(input$groupVar,{
                obj <- objectList$objects[[as.numeric(input$object)]]
                varNames <- extractPhenoData(obj)[[input$groupVar]]
                updateSelectInput(inputId = "blankname", choices = as.character(unique(varNames)))
                updateSelectInput(inputId = "qcname", choices = as.character(unique(varNames)))
                updateSelectInput(inputId = "groupName", choices = as.character(unique(varNames)))
                updateSelectInput(inputId = "filtName", choices = as.character(c("none", unique(varNames))))

            }, ignoreInit = TRUE, ignoreNULL = TRUE)
            observeEvent(input$featProcBut, {
                obj <- objectList$objects[[as.numeric(input$object)]]
                if(input$omic == "Metabolomics"){
                    if(input$prePro){
                        if(any(input$preProFuns == "mvImp")){
                            validate(need(all(apply(extractData(obj),2, function(x){
                                length(which(is.na(x)))/length(x)
                            }) < 0.5) , message = "Missing values from one sample are above 50% use a feature/sample NA filter first"))
                        }
                    }

                    returnData$object <- metabFeatureFilter(features = obj, groupvar = input$groupVar,
                                                        blankfilt = input$blankFilt, blankFoldChange = input$blankRatio,
                                                        blankname = input$blankname, samplename = input$groupName,
                                                        cvqcfilt = input$cvFilt, cvqc_thr = input$maxRDS,
                                                        qcname = input$qcname, nafilter = input$naFilt,
                                                        naratioThr = input$nathr, naratioMethod = input$naMethod,
                                                        varfilter = input$metVarFilt, varfun = input$varfun,
                                                        varthr = input$varthr, varquant = input$quant, intfilter = input$intFilt,
                                                        intensitythr = input$intthr, ism0 = input$ism0Filt, hasan = input$annoFilt,
                                                        sampfilter = input$sampleFilt, maxmv = input$maxNA, filtername = input$filtName,
                                                        prepro = input$prePro, preprofuns = input$preProFuns, mvimpmethod = input$mvImpMeth)
                } else if(input$omic == "Transcriptomics"){
                    returnData$object <- geneFeatureFilter(features = input$obj, entrez = input$entrez, rem.dupEntrez = input$entrez,
                                                    varfilt = input$varfilt, varcutoff = input$quantvar, var.func = input$varfun )
                } else if(input$omic == "Both"){
                    if(compFilter){
                        returnData$object <- groupFeatureFilter(features = input$obj, comptable = input$compTable, pvalthr = input$pthr,
                                                            logFCthr = input$fcthr, padjusted = input$padj)
                    }
                    if(biofilter){
                        if(length(returnData$object) == 0){
                            returnData$object <- featureSelection(features = input$obj, biosigndata = input$bioSign,
                                                            model = input$bioSignMod, scoremin = input$minScore)
                        } else{
                            returnData$object <- featureSelection(features = returnData$object, biosigndata = input$bioSign,
                                                            model = input$bioSignMod, scoremin = input$minScore)
                        }

                    }
                }
                sendSweetAlert(title = "Feature Filter",
                                text = "Your features were filtered successfully!",
                                type = "success", session = session)
                returnData$objectNames <- paste0(as.character(names(objectList$objects[as.numeric(input$object)])), "_f")
                returnData$trigger <- returnData$trigger + 1
            }, ignoreInit = TRUE, ignoreNULL = TRUE)
            return(returnData)
        }
    )
}
