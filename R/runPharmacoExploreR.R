#' Launch PharmacoExploreR Shiny App
#'
#' Launches the interactive Shiny application for exploring pharmacogenomic data.
#' The app allows users to upload PharmacoSet data or use the included demo dataset
#' to perform correlation analysis, visualize results, and classify samples.
#'
#' @return No return value, opens the Shiny app in a web browser.
#'
#' @examples
#' \dontrun{
#' # Launch the Shiny app
#' runPharmacoExploreR()
#' }
#'
#' @export
#' @import shiny
runPharmacoExploreR <- function() {
  app_dir <- system.file("shiny-scripts", package = "PharmacoExploreR")
  
  if (app_dir == "") {
    stop("Could not find Shiny app directory. Try re-installing PharmacoExploreR.")
  }
  
  shiny::runApp(app_dir, display.mode = "normal")
}