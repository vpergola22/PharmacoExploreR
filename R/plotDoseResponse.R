# TODO Documentation
#
#
#
# @param 
#
# @return TODO
# 
#
# @examples TODO
#
#
# @references
#
#
# @export
# @import PharmacoGx
# @import ggplot2
# @import reshape2
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