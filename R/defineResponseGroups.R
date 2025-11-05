#' Define Response Groups Based on Drug Sensitivity
#'
#' Classifies samples (cell lines) into "sensitive" and "resistant" groups
#' based on their drug sensitivity values. Multiple classification methods
#' are supported: median split, quantile-based, or manual threshold.
#'
#' @param pset A PharmacoSet object from the PharmacoGx package.
#' @param drug A character string specifying the drug name for classification.
#' @param sensitivity.measure A character string specifying the sensitivity
#'   metric. Default is "auc_recomputed". For AUC/AAC, lower values typically
#'   indicate greater sensitivity.
#' @param method A character string specifying the classification method.
#'   Options are:
#'   \itemize{
#'     \item "median" - Split at median (50th percentile)
#'     \item "quantile" - Use top and bottom quantiles (default: 25th and 75th)
#'     \item "manual" - Use a user-specified threshold
#'   }
#' @param threshold A numeric value specifying the manual threshold (only used
#'   when method = "manual"). Samples with values <= threshold are classified
#'   as "sensitive".
#' @param quantiles A numeric vector of length 2 specifying the lower and upper
#'   quantiles for classification (only used when method = "quantile").
#'   Default is c(0.25, 0.75). Samples in the middle are classified as NA.
#'
#' @return A named factor vector with levels "sensitive" and "resistant",
#'   where names correspond to sample identifiers. When method = "quantile",
#'   some samples may be NA (middle quantile).
#'
#' @examples
#' \dontrun{
#' library(PharmacoGx)
#' 
#' # Load data
#' pset <- downloadPSet("NCI60_2021")
#' drug_name <- drugNames(pset)[1]
#' 
#' # Method 1: Median split (most common)
#' groups_median <- defineResponseGroups(
#'   pset = pset,
#'   drug = drug_name,
#'   sensitivity.measure = "aac_recomputed",
#'   method = "median"
#' )
#' table(groups_median)
#' 
#' # Method 2: Quantile-based (more stringent)
#' groups_quantile <- defineResponseGroups(
#'   pset = pset,
#'   drug = drug_name,
#'   sensitivity.measure = "aac_recomputed",
#'   method = "quantile",
#'   quantiles = c(0.25, 0.75)
#' )
#' table(groups_quantile, useNA = "ifany")
#' 
#' # Method 3: Manual threshold
#' groups_manual <- defineResponseGroups(
#'   pset = pset,
#'   drug = drug_name,
#'   sensitivity.measure = "aac_recomputed",
#'   method = "manual",
#'   threshold = 0.6
#' )
#' table(groups_manual)
#' }
#'
#' @references
#' Smirnov, P., Safikhani, Z., El-Hachem, N., Wang, D., She, A., Olsen, C.,
#' Freeman, M., Selby, H., Gendoo, D. M., Grossman, P., Beck, A. H.,
#' Aerts, H. J., Lupien, M., Goldenberg, A., & Haibe-Kains, B. (2016).
#' PharmacoGx: an R package for analysis of large pharmacogenomic datasets.
#' \emph{Bioinformatics}, 32(8), 1244-1246.
#' \href{https://doi.org/10.1093/bioinformatics/btv723}{Link}.
#'
#' @export
#' @importFrom PharmacoGx summarizeSensitivityProfiles
#' @importFrom stats median quantile
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