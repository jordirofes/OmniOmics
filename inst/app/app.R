library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(shinyFiles)
library(OmniOmics)
library(xcms)
library(SummarizedExperiment)
library(CAMERA)
library(cliqueMS)
library(pmp)
library(biosigner)
library(genefilter)
library(plotly)
library(ggplot2)
library(caret)
library(biosigner)
library(gplots)
library(DT)
library(limma)
library(ggpubr)
library(oligo)
library(pvca)
library(pd.mogene.2.1.st)
library(ROCR)
library(mogene21sttranscriptcluster.db)
library(RColorBrewer)

i <- list.files(system.file("app", package = "OmniOmics"), full.names = TRUE, pattern = "UI")
i <- c(i, list.files(system.file("app", package = "OmniOmics"), full.names = TRUE, pattern = "Server"))
for(j in 1:length(i)){
    source(i[j], local = TRUE)
}
# source("./inst/importServer.R", )
# source("./inst/importUI.R")
# source("./inst/objectDataServer.R")
# source("./inst/objectDataUI.R")
# source("./inst/processServer.R")
# source("./inst/processUI.R")

ui <- dashboardPage(
    dashboardHeader(title = "OmniOmics"),
    dashboardSidebar(
        sidebarMenu(
            menuItem("Data Import", tabName = "importData", icon = icon("th")),
            menuItem("Spectra Visualization", tabName = "metabVis", icon = icon("th")),
            menuItem("Object Data", tabName = "objectData", icon = icon("th")),
            menuItem("Data Processing", tabName = "procData", icon = icon("th")),
            menuItem("Batch Correction", tabName = "batchCorr", icon = icon("th")),
            menuItem("Feature Processing", tabName = "featureProc", icon = icon("th")),
            menuItem("Group Comparisson", tabName = "groupComp", icon = icon("th")),
            menuItem("Multivariate Analysis", tabName = "multAnal", icon = icon("th")),
            menuItem("Machine Learning", tabName = "machLearn", icon = icon("th"))
    )),
    dashboardBody(
        tabItems(
            tabItem(tabName = "importData", importUI("import")),
            tabItem(tabName = "metabVis", metaboVisUI("metVis")),
            tabItem(tabName = "objectData", objectDataUI("objectDt")),
            tabItem(tabName = "procData", processUI("proc")),
            tabItem(tabName = "batchCorr", batchUI("corr")),
            tabItem(tabName = "featureProc", featureProcUI("featProc")),
            tabItem(tabName = "groupComp", groupCompUI("gComp")),
            tabItem(tabName = "multAnal", multiAnalyUI("mAnal")),
            tabItem(tabName = "machLearn", machineLearnUI("ml"))
        )
    )
)


server <- function(input, output, session){

    objectList <- reactiveValues(objects = list())

    returnImport <- importServer("import", objectList)
    observeEvent(returnImport$trigger,{
        if(class(returnImport$object) == "list"){
            objectListed <- isolate(returnImport$object)
        } else{
            objectListed <- list(isolate(returnImport$object))
        }
        names(objectListed) <- isolate(returnImport$objectName)
        objectList$objects <- c(isolate(objectList$objects), objectListed)
        # objectList$len <- length(isolate(objectList$objects))
    }, ignoreInit = TRUE, ignoreNULL = TRUE, priority = 1)

    objectDataServer("objectDt", objectList)
    metaboVisServer("metVis", objectList)

    returnProc <- processServer("proc", objectList)
    observeEvent(returnProc$trigger, {

        if(class(returnProc$object) == "list"){
            objectListed <- isolate(returnProc$object)
        } else{
            objectListed <- list(isolate(returnProc$object))
        }
        names(objectListed) <- isolate(returnProc$objectNames)
        objectList$objects <- c(isolate(objectList$objects), objectListed)
    }, ignoreInit = TRUE, ignoreNULL = TRUE, priority = 1)

    returnBatch <- batchServer("corr", objectList)
    observeEvent(returnBatch$trigger, {

        objectListed <- list(isolate(returnBatch$object))
        names(objectListed) <- isolate(returnBatch$objectNames)
        objectList$objects <- c(isolate(objectList$objects), objectListed)
    }, ignoreInit = TRUE, ignoreNULL = TRUE, priority = 1)

    returnFeatureProc <- featureProcServer("featProc", objectList)
    observeEvent(returnFeatureProc$trigger, {
        objectListed <- list(isolate(returnFeatureProc$object))
        names(objectListed) <- isolate(returnFeatureProc$objectNames)
        objectList$objects <- c(isolate(objectList$objects), objectListed)
    }, ignoreInit = TRUE, ignoreNULL = TRUE, priority = 1)

    returnGroupComp <- groupCompServer("gComp", objectList)
    observeEvent(returnGroupComp$trigger, {
        if(class(returnGroupComp$object) == "list"){
            objectListed <- isolate(returnGroupComp$object)
        } else{
            objectListed <- list(isolate(returnGroupComp$object))
        }
        names(objectListed) <- isolate(returnGroupComp$objectNames)
        objectList$objects <- c(isolate(objectList$objects), objectListed)
    }, ignoreInit = TRUE, ignoreNULL = TRUE, priority = 1)

    returnMultiAnaly <- multiAnalyServer("mAnal", objectList)
    observeEvent(returnMultiAnaly$trigger, {
        objectListed <- list(isolate(returnMultiAnaly$object))
        names(objectListed) <- isolate(returnMultiAnaly$objectNames)
        objectList$objects <- c(isolate(objectList$objects), objectListed)
    }, ignoreInit = TRUE, ignoreNULL = TRUE, priority = 1)

    returnMachineLearn <- machineLearnServer("ml", objectList)
    observeEvent(returnMachineLearn$trigger, {

        if(class(returnMachineLearn$object) == "list"){
            objectListed <- isolate(returnMachineLearn$object)
        } else{
            objectListed <- list(isolate(returnMachineLearn$object))
        }
        names(objectListed) <- isolate(returnMachineLearn$objectNames)
        objectList$objects <- c(isolate(objectList$objects), objectListed)
    }, ignoreInit = TRUE, ignoreNULL = TRUE, priority = 1)
}


shinyApp(ui, server)
