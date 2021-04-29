library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(shinyFiles)
library(OmniOmics)
library(xcms)

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
    procVars <- processServer("proc")

    objectDataServer("objectDt", objectList)
    print(objectList)
}


shinyApp(ui, server)
