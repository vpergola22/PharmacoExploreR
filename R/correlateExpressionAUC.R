#' Correlate Gene Expression with Drug Response
#'
#' Computes correlations between gene expression levels and drug sensitivity
#' metrics (e.g., AUC, AAC, IC50) for a specified drug across samples in a
#' PharmacoSet object. For each gene, the function calculates the correlation
#' coefficient, p-value, and FDR-adjusted p-value to identify genes whose
#' expression is associated with drug response.
#'
#' @param pset A PharmacoSet object from the PharmacoGx package containing
#'   molecular profiling and drug sensitivity data.
#' @param drug A character string specifying the drug name to analyze. Must
#'   match a drug name in the sensitivity data. Use \code{drugNames(pset)}
#'   to see available options.
#' @param mDataType A character string specifying the molecular data type.
#'   Default is "rna" for gene expression. Other options depend on the
#'   PharmacoSet (e.g., "mutation", "cnv").
#' @param sensitivity.measure A character string specifying the sensitivity
#'   metric to use. Default is "auc_recomputed". Common options include
#'   "aac_recomputed", "ic50_published", "auc_published".
#' @param method A character string specifying the correlation method.
#'   Options are "pearson" (default), "spearman", or "kendall".
#'
#' @return A data frame with the following columns:
#' \itemize{
#'   \item gene - Character vector of gene names
#'   \item cor - Numeric vector of correlation coefficients
#'   \item pval - Numeric vector of p-values from correlation tests
#'   \item adj_pval - Numeric vector of FDR-adjusted p-values
#' }
#'
#' @examples
#' # Load the bundled mini PharmacoSet dataset
#' data("nci60_mini")
#'
#' # View available drugs
#' head(PharmacoGx::drugNames(nci60_mini))
#'
#' # Correlate gene expression with drug response
#' results <- correlateExpressionAUC(
#'   pset = nci60_mini,
#'   drug = "Erlotinib",
#'   mDataType = "rna",
#'   sensitivity.measure = "aac_recomputed",
#'   method = "pearson"
#' )
#'
#' # View top correlated genes
#' head(results[order(abs(results$cor), decreasing = TRUE), ])
#'
#' # Filter for significant genes (safe even if adj_pval has NAs)
#' sig_genes <- subset(results, !is.na(adj_pval) && adj_pval < 0.05)
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
#' @export
#' @importFrom PharmacoGx summarizeMolecularProfiles summarizeSensitivityProfiles
#' @importFrom SummarizedExperiment assay
#' @importFrom stats cor.test p.adjust complete.cases var
#' @importFrom utils head
correlateExpressionAUC <- function(pset, drug, mDataType = "rna", 
                                   sensitivity.measure = "auc_recomputed", 
                                   method = "pearson"){
  # ensure PharmacoGx is available
  if (!requireNamespace("PharmacoGx", quietly = TRUE)) {
    stop("PharmacoGx package required.")
  }
  
  # check for valid PSet
  if (!inherits(pset, "PharmacoSet")) {
    stop("Input must be a PharmacoSet object.")
  }
  
  # extract expression and sensitivity data
  expr <- PharmacoGx::summarizeMolecularProfiles(pset, mDataType = mDataType, 
                                                 summary.stat = "mean")
  auc <- PharmacoGx::summarizeSensitivityProfiles(pset, 
                                                  sensitivity.measure = sensitivity.measure)
  
  # extract expression matrix from SummarizedExperiment if needed
  if (inherits(expr, "SummarizedExperiment")) {
    expr <- SummarizedExperiment::assay(expr)
  }
  
  # get available drugs
  available_drugs <- rownames(auc)
  
  # check that drug is found in the sensitivity data
  if (!(drug %in% available_drugs)) {
    message("\nThe drug '", drug, "' was not found in this PSet.")
    message("Here are some available drugs in this dataset:\n")
    print(utils::head(available_drugs, 15))
    stop("\nPlease select one of the above drug names.")
  }
  
  # ensure sample alignment
  common_samples <- intersect(colnames(expr), colnames(auc))
  if (length(common_samples) == 0) {
    stop("No overlapping samples between expression and sensitivity data.")
  }
  
  # subset to common samples
  expr <- expr[, common_samples, drop = FALSE]
  auc_vals <- auc[drug, common_samples]
  
  # ensure numeric
  auc_vals <- as.numeric(auc_vals)
  
  # check for valid data
  if (all(is.na(auc_vals))) {
    stop("All AUC values are NA for drug: ", drug)
  }
  
  # initialize results matrix
  n_genes <- nrow(expr)
  cor_res <- matrix(NA_real_, nrow = n_genes, ncol = 2,
                    dimnames = list(rownames(expr), c("cor", "pval")))
  
  # compute correlations
  for (i in seq_len(n_genes)) {
    expr_vals <- as.numeric(expr[i, ])
    res <- .safeCor(expr_vals, auc_vals, method)
    cor_res[i, ] <- res
  }
  
  # convert to data frame
  cor_df <- as.data.frame(cor_res, stringsAsFactors = FALSE)
  cor_df$gene <- rownames(cor_res)
  cor_df$adj_pval <- p.adjust(cor_df$pval, method = "fdr")
  
  # reorder columns
  cor_df <- cor_df[, c("gene", "cor", "pval", "adj_pval")]
  rownames(cor_df) <- NULL
  
  return(cor_df)
}






#' Helper Function for Safe Correlation Computation
#'
#' Internal function that safely computes correlation between two numeric
#' vectors, handling missing values and edge cases. Not intended for direct
#' user access.
#'
#' @param x A numeric vector.
#' @param y A numeric vector of the same length as x.
#' @param method Correlation method: "pearson", "spearman", or "kendall".
#'
#' @return A named numeric vector with two elements: "cor" (correlation
#'   coefficient) and "pval" (p-value). Returns NA for both if computation
#'   fails.
#'
#' @keywords internal
#' @noRd
.safeCor <- function(x, y, method = "pearson") {
  # Ensure both numeric
  if (!is.numeric(x) || !is.numeric(y)) {
    return(c(cor = NA_real_, pval = NA_real_))
  }
  
  # remove missing values
  ok <- stats::complete.cases(x, y)
  if (sum(ok) < 3) {
    return(c(cor = NA_real_, pval = NA_real_))
  }
  
  x <- x[ok]
  y <- y[ok]
  
  # handle zero variance
  if (stats::var(x) == 0 || stats::var(y) == 0) {
    return(c(cor = NA_real_, pval = NA_real_))
  }
  
  # compute correlation
  res <- tryCatch({
    test <- suppressWarnings(stats::cor.test(x, y, method = method))
    c(cor = unname(test$estimate), pval = unname(test$p.value))
  }, error = function(e) {
    c(cor = NA_real_, pval = NA_real_)
  })
  
  return(res)
}