#' Plot Gene Expression by Response Group
#'
#' Creates a boxplot or violin plot showing the distribution of a specific
#' gene's expression across drug-sensitive and drug-resistant sample groups.
#' Useful for visualizing differential expression results and validating
#' biomarker candidates.
#'
#' @param pset A PharmacoSet object from the PharmacoGx package.
#' @param gene A character string specifying the gene name to plot.
#' @param groupLabels A named factor vector with levels "sensitive" and
#'   "resistant", typically output from \code{defineResponseGroups()}.
#' @param mDataType A character string specifying the molecular data type.
#'   Default is "rna" for gene expression.
#' @param plotType A character string specifying the plot type. Options are
#'   "boxplot" (default) or "violin".
#' @param showPoints A logical indicating whether to overlay individual
#'   data points. Default is TRUE.
#' @param pointColor A character string specifying the color for individual
#'   points. Default is "#0072B2".
#'
#' @return A ggplot2 object displaying the gene expression distribution
#'   by response group.
#'
#' @examples
#' \dontrun{
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
#' # Run differential expression
#' diff_expr <- runDiffExpr(pset, groups, method = "t.test")
#' 
#' # Plot top gene with boxplot
#' top_gene <- diff_expr$gene[which.min(diff_expr$adj_pval)]
#' 
#' p1 <- plotGeneBoxplot(
#'   pset = nci60_mini,
#'   gene = top_gene,
#'   groupLabels = groups,
#'   plotType = "boxplot"
#' )
#' print(p1)
#' 
#' # Plot with violin plot
#' p2 <- plotGeneBoxplot(
#'   pset = nci60_mini,
#'   gene = top_gene,
#'   groupLabels = groups,
#'   plotType = "violin",
#'   showPoints = TRUE
#' )
#' print(p2)
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
#' @importFrom PharmacoGx summarizeMolecularProfiles
#' @importFrom SummarizedExperiment assay
runDiffExpr <- function(pset, groupLabels, mDataType = "rna", method = c("t.test", "limma")) {
  
  # check inputs
  if (!inherits(pset, "PharmacoSet")) stop("Input must be a PharmacoSet object.")
  if (!is.factor(groupLabels) || !all(c("sensitive", "resistant") %in% levels(groupLabels))) {
    stop("groupLabels must be a factor with levels 'sensitive' and 'resistant'.")
  }
  
  method <- match.arg(method)
  
  # extract expression data
  expr <- PharmacoGx::summarizeMolecularProfiles(pset, mDataType = mDataType, summary.stat = "mean")
  if (inherits(expr, "SummarizedExperiment")) expr <- SummarizedExperiment::assay(expr)
  expr <- as.matrix(expr)
  
  # align samples between expression data and group labels
  common_samples <- intersect(colnames(expr), names(groupLabels))
  if (length(common_samples) == 0) stop("No overlapping samples between expression data and group labels.")
  
  expr <- expr[, common_samples, drop = FALSE]
  groups <- droplevels(groupLabels[common_samples])
  
  # initialize results data frame
  res <- data.frame(
    gene = rownames(expr),
    logFC = NA_real_,
    pval = NA_real_,
    adj_pval = NA_real_,
    stringsAsFactors = FALSE
  )
  
  if (method == "t.test") {
    # compute t-tests for each gene
    for (i in seq_len(nrow(expr))) {
      x <- expr[i, groups == "sensitive"]
      y <- expr[i, groups == "resistant"]
      if (length(x) < 2 || length(y) < 2) {
        res$logFC[i] <- NA
        res$pval[i] <- NA
      } else {
        res$logFC[i] <- mean(x) - mean(y)
        t_res <- try(t.test(x, y), silent = TRUE)
        if (inherits(t_res, "try-error")) {
          res$pval[i] <- NA
        } else {
          res$pval[i] <- t_res$p.value
        }
      }
    }
  } else if (method == "limma") {
    # use limma for differential expression
    if (!requireNamespace("limma", quietly = TRUE)) stop("Package 'limma' is required for method = 'limma'")
    design <- model.matrix(~ groups)
    fit <- limma::lmFit(expr, design)
    fit <- limma::eBayes(fit)
    top <- limma::topTable(fit, coef = 2, number = nrow(expr), sort.by = "none")
    res$logFC <- top$logFC
    res$pval <- top$P.Value
  }
  
  # adjust p-values for multiple testing
  res$adj_pval <- p.adjust(res$pval, method = "fdr")
  
  return(res)
}