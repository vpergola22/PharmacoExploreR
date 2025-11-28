#' Create Scatterplot of Gene Expression vs Drug Response
#'
#' Generates a scatterplot showing the relationship between a specific gene's
#' expression and drug sensitivity (AUC) across samples. The plot includes a
#' linear regression line and displays the correlation coefficient and
#' statistical significance.
#'
#' @param pset A PharmacoSet object from the PharmacoGx package.
#' @param corResults A data frame output from \code{correlateExpressionAUC()},
#'   containing correlation results for genes.
#' @param gene A character string specifying the gene name to plot.
#' @param drug A character string specifying the drug name to plot.
#' @param method A character string specifying the correlation method used.
#'   Options are "pearson" (default), "spearman", or "kendall".
#' @param mDataType A character string specifying the molecular data type.
#'   Default is "rna".
#' @param sensitivity.measure A character string specifying the sensitivity
#'   metric. Default is "auc_recomputed".
#'
#' @return A ggplot2 object displaying the scatterplot with regression line.
#'
#' @examples
#' \dontrun{
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
#' # Plot top correlated gene
#' top_gene <- results$gene[which.max(abs(results$cor))]
#' 
#' plot <- plotExprAUC(
#'   pset = nci60_mini,
#'   corResults = results,
#'   gene = top_gene,
#'   drug = drug_name,
#'   sensitivity.measure = "aac_recomputed"
#' )
#' 
#' print(plot)
#' }
#'
#' @references
#' Morgan, M., Obenchain, V., Hester, J., Pagès, H. (2024).
#' SummarizedExperiment: A container (S4 class) for matrix-like assays.
#' R package version 1.36.0.
#' \href{https://doi.org/10.18129/B9.bioc.SummarizedExperiment}{Link}.
#' 
#' Smirnov, P., Safikhani, Z., El-Hachem, N., Wang, D., She, A., Olsen, C.,
#' Freeman, M., Selby, H., Gendoo, D. M., Grossman, P., Beck, A. H.,
#' Aerts, H. J., Lupien, M., Goldenberg, A., & Haibe-Kains, B. (2016).
#' PharmacoGx: an R package for analysis of large pharmacogenomic datasets.
#' \emph{Bioinformatics}, 32(8), 1244-1246.
#' \href{https://doi.org/10.1093/bioinformatics/btv723}{Link}.
#' 
#' Wickham, H. (2016). ggplot2: Elegant Graphics for Data Analysis.
#' Springer-Verlag New York. \href{https://ggplot2.tidyverse.org}{Link}.
#'
#' @export
#' @import ggplot2
#' @importFrom PharmacoGx summarizeMolecularProfiles summarizeSensitivityProfiles
#' @importFrom SummarizedExperiment assay
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