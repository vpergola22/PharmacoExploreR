#' Plot Dose-Response Curves for Cell Lines
#'
#' Generates dose-response curves showing cell viability (or other response
#' measure) as a function of drug concentration for specified cell lines.
#' Each curve represents one cell line, allowing visual comparison of drug
#' sensitivity across samples.
#'
#' @param pset A PharmacoSet object from the PharmacoGx package containing
#'   raw dose-response data.
#' @param cell.lines A character vector specifying which cell lines to plot.
#'   If NULL, plots all available cell lines (not recommended for large
#'   datasets). Use \code{cellNames(pset)} to see available options.
#' @param sensitivity.measure A character string specifying which response
#'   measure to plot. Default is "Viability". Other options depend on the
#'   PSet and may include measures like "Dose".
#'
#' @return A ggplot2 object displaying dose-response curves with log-scaled
#'   x-axis (dose) and smoothed trend lines for each cell line.
#'
#' @examples
#' library(PharmacoGx)
#' 
#' # Load data
#' # pset <- downloadPSet("NCI60_2021")
#' 
#' # (for demo purposes we will use the mini dataset)
#' data("nci60_mini")
#' 
#' # Extract the cell lines that actually have dose-response data
#' valid_lines <- rownames(nci60_mini@@treatmentResponse$raw[, , "Viability"])
#' selected_lines <- valid_lines[1:3]
#'
#' # Plot dose-response curves
#' plot <- plotDoseResponse(
#'   pset = nci60_mini,
#'   cell.lines = selected_lines
#' )
#' print(plot)
#'
#' @references
#' Smirnov, P., Safikhani, Z., El-Hachem, N., Wang, D., She, A., Olsen, C.,
#' Freeman, M., Selby, H., Gendoo, D. M., Grossman, P., Beck, A. H.,
#' Aerts, H. J., Lupien, M., Goldenberg, A., & Haibe-Kains, B. (2016).
#' PharmacoGx: an R package for analysis of large pharmacogenomic datasets.
#' \emph{Bioinformatics}, 32(8), 1244-1246.
#' \href{https://doi.org/10.1093/bioinformatics/btv723}{Link}.
#' 
#' Wickham, H. (2007). Reshaping Data with the reshape Package.
#' \emph{Journal of Statistical Software}, 21(12), 1-20.
#' \href{https://www.jstatsoft.org/v21/i12/}{Link}.
#'
#' Wickham, H. (2016). ggplot2: Elegant Graphics for Data Analysis.
#' Springer-Verlag New York. \href{https://ggplot2.tidyverse.org}{Link}.
#'
#' @export
#' @import ggplot2
#' @importFrom reshape2 melt
plotDoseResponse <- function(pset, cell.lines = NULL, sensitivity.measure = "Viability") {
  
  # check valid sensitivity.measure
  valid_measures <- dimnames(pset@treatmentResponse$raw)[[3]]
  if (!(sensitivity.measure %in% valid_measures)) {
    stop("Invalid sensitivity.measure. Choose among: ", paste(valid_measures, collapse = ", "))
  }
  
  # extract dose and response matrices
  dose_mat <- pset@treatmentResponse$raw[, , "Dose"]
  resp_mat <- pset@treatmentResponse$raw[, , sensitivity.measure]
  
  # if user specified cell lines, filter
  if (!is.null(cell.lines)) {
    missing_lines <- setdiff(cell.lines, rownames(resp_mat))
    if (length(missing_lines) > 0) {
      stop("The following cell lines are not in this PSet: ", paste(missing_lines, collapse = ", "))
    }
    dose_mat <- dose_mat[cell.lines, , drop = FALSE]
    resp_mat <- resp_mat[cell.lines, , drop = FALSE]
  }
  
  # long format
  long_df <- melt(resp_mat, varnames = c("CellLine", "DoseIdx"), value.name = "Response")
  long_df$Dose <- as.numeric(melt(dose_mat, varnames = c("CellLine", "DoseIdx"))$value)
  
  # remove NAs
  long_df <- na.omit(long_df)
  
  if (nrow(long_df) == 0) stop("No valid data available for plotting.")
  
  # plot
  ggplot(long_df, aes(x = Dose, y = Response, color = CellLine)) +
    geom_point() +
    geom_smooth(method = "loess", se = FALSE) +
    scale_x_log10() +
    labs(x = "Dose", y = sensitivity.measure, title = "Dose-Response Curves") +
    theme_minimal() +
    theme(legend.position = "bottom")
}