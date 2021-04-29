library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(shinyFiles)
library(OmniOmics)

ui <- dashboardPage(
    dashboardHeader(title = "OmniOmics"),
    dashboardSidebar(
        sidebarMenu(
            menuItem("Data Import", tabName = "importData", icon = icon("th")),
            menuItem("Data Processing", tabName = "procData", icon = icon("th"))
    )),
    dashboardBody(
        tabItems(
            tabItem(tabName = "importData",
                importUI("import")
            ),
            tabItem(tabName = "procData",
                processUI("proc")
            )
        )
    )
)


server <- function(input, output){
    importVars <- importServer("import")
    procVars <- processServer("proc")


}


shinyApp(ui, server)
