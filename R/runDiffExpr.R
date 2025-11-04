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