#' code for running shiny app
#'
#' @import shiny shinydashboard shinyjs DT plotly shinyBS
#' @export
run_nMyo <- function() {
  require(shinydashboard)
  require(shiny)
  require(shinyjs)
  require(DT)
  require(shinyBS)
  require(plotly)
  appDir <- system.file("shiny-examples", "nmyo_app", package = "nMyo")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `nMyo`.", call. = FALSE)
  }

  runApp(appDir, display.mode = "normal")
}
