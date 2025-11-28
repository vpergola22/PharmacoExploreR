#' Create Volcano Plot of Gene-Drug Correlations
#'
#' Generates a volcano plot displaying the correlation strength (x-axis)
#' versus statistical significance (y-axis) for all genes. This provides
#' a genome-wide overview of gene-drug associations, highlighting genes
#' that meet both correlation and significance thresholds.
#'
#' @param corResults A data frame output from \code{correlateExpressionAUC()},
#'   containing columns: gene, cor (correlation), pval, and adj_pval.
#' @param corThreshold A numeric value specifying the absolute correlation
#'   threshold for significance. Default is 0.3.
#' @param pvalThreshold A numeric value specifying the adjusted p-value
#'   threshold for significance. Default is 0.05.
#'
#' @return A ggplot2 object displaying the volcano plot with threshold lines
#'   and colored points indicating significant genes.
#'
#' @examples
#' library(PharmacoGx)
#' 
#' # Load data and compute correlations
#' # pset <- downloadPSet("NCI60_2021")
#' 
#' # (for demo purposes we will use the mini dataset)
#' data("nci60_mini")
#' 
#' drug_name <- drugNames(nci60_mini)[1]
#' 
#' results <- correlateExpressionAUC(
#'   pset = nci60_mini,
#'   drug = drug_name,
#'   sensitivity.measure = "aac_recomputed"
#' )
#' 
#' # Create volcano plot
#' volcano <- volcanoAUC(
#'   corResults = results,
#'   corThreshold = 0.3,
#'   pvalThreshold = 0.05
#' )
#' 
#' print(volcano)
#' 
#' # Use stricter thresholds
#' volcano_strict <- volcanoAUC(
#'   corResults = results,
#'   corThreshold = 0.5,
#'   pvalThreshold = 0.01
#' )
#'
#' @references
#' Wickham, H. (2016). ggplot2: Elegant Graphics for Data Analysis.
#' Springer-Verlag New York. \href{https://ggplot2.tidyverse.org}{Link}.
#'
#' @export
#' @import ggplot2

volcanoAUC <- function(corResults, corThreshold = 0.3, pvalThreshold = 0.05) {
  # ensure expected columns
  if (!all(c("gene", "cor", "adj_pval") %in% colnames(corResults))) {
    stop("Input data frame must contain 'gene', 'cor', and 'adj_pval' columns.")
  }
  
  df <- corResults
  df$significant <- with(df, abs(cor) >= corThreshold & adj_pval <= pvalThreshold)
  
  ggplot(df, aes(x = cor, y = -log10(adj_pval))) +
    geom_point(aes(color = significant), alpha = 0.7) +
    scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "#D55E00")) +
    geom_vline(xintercept = c(-corThreshold, corThreshold), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(pvalThreshold), linetype = "dashed", color = "black") +
    theme_minimal(base_size = 13) +
    labs(
      title = "Volcano Plot of Gene–AUC Correlations",
      x = "Correlation Coefficient (r)",
      y = expression(-log[10](Adjusted~p~value)),
      color = "Significant"
    ) +
    theme(plot.title = element_text(face = "bold"))
}