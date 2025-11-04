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
plotExprAUC <- function(pset, corResults, gene, drug, method = "pearson", mDataType = "rna", sensitivity.measure = "auc_recomputed") {
  
  if (!inherits(pset, "PharmacoSet")) {
    stop("Input must be a PharmacoSet object.")
  }
  
  if (missing(gene) || missing(drug)) {
    stop("Please specify both 'gene' and 'drug'.")
  }
  
  # Extract molecular and sensitivity data
  expr <- PharmacoGx::summarizeMolecularProfiles(pset, mDataType = mDataType, summary.stat = "mean")
  auc <- PharmacoGx::summarizeSensitivityProfiles(pset, sensitivity.measure = sensitivity.measure)
  
  # If SummarizedExperiment, extract assay matrix
  if (inherits(expr, "SummarizedExperiment")) {
    expr <- SummarizedExperiment::assay(expr)
  }
  expr <- as.matrix(expr)
  
  # Extract numeric AUC vector robustly
  if (drug %in% rownames(auc)) {
    auc_vals <- as.numeric(auc[drug, , drop = TRUE])
    names(auc_vals) <- colnames(auc)
  } else if (drug %in% colnames(auc)) {
    auc_vals <- as.numeric(auc[, drug, drop = TRUE])
    names(auc_vals) <- rownames(auc)
  } else {
    stop(paste("Drug", drug, "not found in AUC data."))
  }
  
  # Align samples exactly as in correlateExpressionAUC()
  common_samples <- intersect(colnames(expr), names(auc_vals))
  if (length(common_samples) == 0) stop("No overlapping samples between expression and AUC data.")
  
  expr_vec <- as.numeric(expr[gene, common_samples])
  auc_vec <- as.numeric(auc_vals[common_samples])
  
  # Prepare data frame
  df <- data.frame(Expression = expr_vec, AUC = auc_vec)
  
  # Retrieve correlation from corResults if available
  cor_row <- corResults[corResults$gene == gene, ]
  if (nrow(cor_row) == 1) {
    cor_label <- sprintf("r = %.2f (adj p = %.2g)", cor_row$cor, cor_row$adj_pval)
  } else {
    # fallback if gene not found
    cor_res <- suppressWarnings(cor.test(df$Expression, df$AUC, method = method))
    cor_label <- sprintf("%s = %.2f (p = %.2g)",
                         ifelse(method == "pearson", "r", "rho"),
                         cor_res$estimate, cor_res$p.value)
  }
  
  # Plot
  ggplot(df, aes(x = Expression, y = AUC)) +
    geom_point(alpha = 0.7, color = "#0072B2") +
    geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
    theme_minimal(base_size = 13) +
    labs(
      title = paste0("Expression of ", gene, " vs AUC for ", drug),
      subtitle = cor_label,
      x = "Gene Expression (normalized)",
      y = "Drug Response (AUC)"
    ) +
    theme(plot.title = element_text(face = "bold"))
}