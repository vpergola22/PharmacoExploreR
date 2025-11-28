#' Perform Differential Expression Analysis Between Response Groups
#'
#' Identifies genes that are differentially expressed between drug-sensitive
#' and drug-resistant samples. Supports t-test or limma-based analysis methods.
#' Returns statistical measures including log fold change and adjusted p-values
#' for each gene.
#'
#' @param pset A PharmacoSet object from the PharmacoGx package.
#' @param groupLabels A named factor vector with levels "sensitive" and
#'   "resistant", typically output from \code{defineResponseGroups()}.
#'   Names must correspond to sample identifiers in the PSet.
#' @param mDataType A character string specifying the molecular data type.
#'   Default is "rna" for gene expression.
#' @param method A character string specifying the statistical method.
#'   Options are:
#'   \itemize{
#'     \item "t.test" - Welch's t-test for each gene (default)
#'     \item "limma" - Linear modeling with empirical Bayes (requires limma package)
#'   }
#'
#' @return A data frame with the following columns:
#' \itemize{
#'   \item gene - Character vector of gene names
#'   \item logFC - Numeric vector of log fold changes (sensitive - resistant)
#'   \item pval - Numeric vector of p-values
#'   \item adj_pval - Numeric vector of FDR-adjusted p-values
#' }
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
#' drug_name <- drugNames(nci60_mini)[1]
#' 
#' # Define response groups
#' groups <- defineResponseGroups(
#'   pset = nci60_mini,
#'   drug = drug_name,
#'   sensitivity.measure = "aac_recomputed",
#'   method = "median"
#' )
#' 
#' # Run differential expression using t-test
#' diff_expr <- runDiffExpr(
#'   pset = nci60_mini,
#'   groupLabels = groups,
#'   mDataType = "rna",
#'   method = "t.test"
#' )
#' 
#' # View top differentially expressed genes
#' head(diff_expr[order(diff_expr$adj_pval), ])
#' 
#' # Filter for significant genes
#' sig_genes <- diff_expr[diff_expr$adj_pval < 0.05 & abs(diff_expr$logFC) > 1, ]
#' 
#' # Using limma (if installed)
#' diff_expr_limma <- runDiffExpr(
#'   pset = nci60_mini,
#'   groupLabels = groups,
#'   method = "limma"
#' )
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
#' R Core Team (2024). R: A Language and Environment for Statistical Computing.
#' R Foundation for Statistical Computing, Vienna, Austria.
#' \href{https://www.R-project.org/}{Link}.
#' 
#' Ritchie, M. E., Phipson, B., Wu, D., Hu, Y., Law, C. W., Shi, W., &
#' Smyth, G. K. (2015). limma powers differential expression analyses for
#' RNA-sequencing and microarray studies. \emph{Nucleic Acids Research},
#' 43(7), e47. \href{https://doi.org/10.1093/nar/gkv007}{Link}.
#'
#' @export
#' @importFrom PharmacoGx summarizeMolecularProfiles
#' @importFrom SummarizedExperiment assay
#' @importFrom stats t.test p.adjust model.matrix
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