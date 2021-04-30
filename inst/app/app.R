library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(shinyFiles)
library(OmniOmics)
library(xcms)
library(SummarizedExperiment)
library(CAMERA)
library(cliqueMS)

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
            menuItem("Object Data", tabName = "objectData", icon = icon("th")),
            menuItem("Data Processing", tabName = "procData", icon = icon("th"))
    )),
    dashboardBody(
        tabItems(
            tabItem(tabName = "importData",
                importUI("import")
            ),
            tabItem(tabName = "objectData",
                objectDataUI("objectDt")
            ),
            tabItem(tabName = "procData",
                processUI("proc")
            )
        )
    )
)


server <- function(input, output, session){

    objectList <- reactiveValues(objects = list())

    returnImport <- importServer("import")
    observeEvent(returnImport$trigger,{
        objectListed <- list(isolate(returnImport$object))
        names(objectListed) <- isolate(returnImport$objectName)
        objectList$objects <- c(isolate(objectList$objects), objectListed)
        # objectList$len <- length(isolate(objectList$objects))
    }, ignoreInit = TRUE, ignoreNULL = TRUE, priority = 1)
    returnProc <- processServer("proc", objectList)
    observeEvent(returnProc$trigger, {
        objectListed <- isolate(returnProc$object)
        names(objectListed) <- isolate(returnProc$objectNames)
        objectList$objects <- c(isolate(objectList$objects), objectListed)
    }, ignoreInit = TRUE, ignoreNULL = TRUE, priority = 1)

    objectDataServer("objectDt", objectList)

    print(objectList)
}


shinyApp(ui, server)
