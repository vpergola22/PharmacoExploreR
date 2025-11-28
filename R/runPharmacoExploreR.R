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
#' @references
#' Chang, W., Cheng, J., Allaire, J., Sievert, C., Schloerke, B., Xie, Y.,
#' Allen, J., McPherson, J., Dipert, A., Borges, B. (2023). shiny: Web
#' Application Framework for R. R package version 1.8.0.
#' \href{https://CRAN.R-project.org/package=shiny}{Link}.
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