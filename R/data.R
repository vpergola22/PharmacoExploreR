#' Mini NCI60 PharmacoSet
#'
#' A subset of the NCI60 pharmacogenomic dataset created for PharmacoExploreR
#' examples and vignettes. Contains realistic gene expression and drug sensitivity
#' data with proper structure to eliminate warnings during analysis.
#' 
#' @source Shoemaker RH. The NCI60 human tumour cell line anticancer drug screen.
#' Nat Rev Cancer. 2006 Oct;6(10):813-23. doi: 10.1038/nrc1951.
#'
#' @docType data
#' @usage data(nci60_mini)
#'
#' @format A PharmacoSet object with:
#' \describe{
#'   \item{molecularProfiles}{Expression data for 150 cell lines}
#'   \item{treatmentResponse}{
#'     Drug sensitivity for 10 drugs including:
#'     \itemize{
#'       \item AAC (Area Above Curve)
#'       \item AUC (Area Under Curve)
#'       \item IC50 values
#'       \item Raw dose-response curves with 8 dose points
#'     }
#'   }
#'   \item{sample}{Metadata for 150 cell lines from 5 tissue types}
#'   \item{treatment}{Metadata for 10 drugs}
#' }
#'
#' @examples
#' \dontrun{
#' data(nci60_mini)
#' library(PharmacoGx)
#' drugNames(nci60_mini)
#' cellNames(nci60_mini)
#' results <- correlateExpressionAUC(nci60_mini, drug = "Erlotinib")
#' head(results)
#' }
"nci60_mini"
