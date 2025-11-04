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






# Helper function TODO documentation
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