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
plotGeneBoxplot <- function(pset, gene, groupLabels, mDataType = "rna",
                            plotType = c("boxplot", "violin"), 
                            showPoints = TRUE, pointColor = "#0072B2") {
  
  # Check inputs
  if (!inherits(pset, "PharmacoSet")) stop("Input must be a PharmacoSet object.")
  if (!gene %in% rownames(PharmacoGx::summarizeMolecularProfiles(pset, mDataType = mDataType))) {
    stop(paste("Gene", gene, "not found in expression data."))
  }
  if (!is.factor(groupLabels) || !all(c("sensitive", "resistant") %in% levels(groupLabels))) {
    stop("groupLabels must be a factor with levels 'sensitive' and 'resistant'.")
  }
  
  plotType <- match.arg(plotType)
  
  # Extract expression data
  expr <- PharmacoGx::summarizeMolecularProfiles(pset, mDataType = mDataType, summary.stat = "mean")
  if (inherits(expr, "SummarizedExperiment")) expr <- SummarizedExperiment::assay(expr)
  expr <- as.matrix(expr)
  
  # Align samples
  common_samples <- intersect(colnames(expr), names(groupLabels))
  expr_vec <- expr[gene, common_samples]
  groups <- droplevels(groupLabels[common_samples])
  
  # Create data frame for plotting
  df <- data.frame(
    Expression = expr_vec,
    Group = groups
  )
  
  # Create plot
  p <- ggplot(df, aes(x = Group, y = Expression, fill = Group)) +
    theme_minimal(base_size = 13) +
    labs(title = paste("Expression of", gene, "by Drug Response"),
         x = "Response Group", y = "Expression") +
    scale_fill_manual(values = c("sensitive" = "#009E73", "resistant" = "#D55E00")) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "none")
  
  if (plotType == "boxplot") {
    p <- p + geom_boxplot(alpha = 0.7)
    if (showPoints) p <- p + geom_jitter(width = 0.2, alpha = 0.6, color = pointColor)
  } else if (plotType == "violin") {
    p <- p + geom_violin(alpha = 0.7)
    if (showPoints) p <- p + geom_jitter(width = 0.2, alpha = 0.6, color = pointColor)
  }
  
  return(p)
}