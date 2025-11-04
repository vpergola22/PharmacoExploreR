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
defineResponseGroups <- function(pset, drug, 
                                 sensitivity.measure = "auc_recomputed", 
                                 method = c("median", "quantile", "manual"),
                                 threshold = NULL,
                                 quantiles = c(0.25, 0.75)) {
  # check inputs
  if (!inherits(pset, "PharmacoSet")) {
    stop("Input must be a PharmacoSet object.")
  }
  
  method <- match.arg(method)
  
  # extract sensitivity data
  auc <- PharmacoGx::summarizeSensitivityProfiles(pset, sensitivity.measure = sensitivity.measure)
  
  # extract numeric vector for the drug
  if (drug %in% rownames(auc)) {
    resp_vals <- as.numeric(auc[drug, , drop = TRUE])
    names(resp_vals) <- colnames(auc)
  } else if (drug %in% colnames(auc)) {
    resp_vals <- as.numeric(auc[, drug, drop = TRUE])
    names(resp_vals) <- rownames(auc)
  } else {
    stop(paste("Drug", drug, "not found in sensitivity data."))
  }
  
  # remove NA values
  resp_vals <- resp_vals[!is.na(resp_vals)]
  
  # determine thresholds
  if (method == "median") {
    thresh <- median(resp_vals)
    group <- ifelse(resp_vals <= thresh, "sensitive", "resistant")
  } else if (method == "quantile") {
    q_low <- quantile(resp_vals, probs = quantiles[1])
    q_high <- quantile(resp_vals, probs = quantiles[2])
    group <- ifelse(resp_vals <= q_low, "sensitive",
                    ifelse(resp_vals >= q_high, "resistant", NA))
  } else if (method == "manual") {
    if (is.null(threshold)) stop("For 'manual' method, provide a numeric threshold.")
    group <- ifelse(resp_vals <= threshold, "sensitive", "resistant")
  }
  
  # return as a named factor vector
  group <- factor(group, levels = c("sensitive", "resistant"))
  return(group)
}