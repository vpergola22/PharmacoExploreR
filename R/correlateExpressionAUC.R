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
correlateExpressionAUC <- function(pset, drug, method = "pearson"){
  # ensure PharmacoGx is available
  if (!requireNamespace("PharmacoGx", quietly = TRUE)) {
    stop("PharmacoGx package required.")
  }
  
  # check for valid PSet
  if (!inherits(pset, "PharmacoSet")) {
    stop("Input must be a PharmacoSet object.")
  }
  
  # Extract expression and sensitivity data
  expr <- PharmacoGx::summarizeMolecularProfiles(pset, mDataType = "rna", summary.stat = "mean") # TODO double chec that these are giving the correct metric
  auc <- PharmacoGx::summarizeSensitivityProfiles(pset, sensitivity.measure = "auc_recomputed")
  
  
  # check that drug is found in the sensitivity data
  if (!(drug %in% colnames(auc))) {
    stop(paste("Drug", drug, "not found in sensitivity data."))
  }
  
  # map samples from expression data to sensitivity data
  common_samples <- intersect(colnames(expr), rownames(auc))
  expr <- expr[, common_samples, drop = FALSE]
  auc_vals <- auc[common_samples, drug]
  
  
  # compute correlation for each gene
  
  
}