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