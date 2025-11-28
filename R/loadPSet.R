#' Load a PharmacoSet from PharmacoGx
#'
#' Downloads and loads a PharmacoSet object from the PharmacoGx repository.
#' PharmacoSets are standardized data structures containing pharmacogenomic
#' data including molecular profiles (e.g., gene expression) and drug
#' sensitivity measurements for cancer cell lines or patient samples.
#'
#' @param psetName A character string specifying the name of the PharmacoSet
#'   to download. Must match an available PSet name from PharmacoGx. Use
#'   \code{PharmacoGx::availablePSets()} to see all options.
#'
#' @return A PharmacoSet object containing molecular profiling data, drug
#'   sensitivity measurements, and associated metadata.
#'
#' @examples
#' # Load the NCI60 dataset
#' pset <- loadPSet("NCI60_2021")
#' 
#' # View available datasets first
#' library(PharmacoGx)
#' availablePSets()
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
#' @importFrom PharmacoGx availablePSets downloadPSet
loadPSet <- function(psetName){
  # ensure PharmacoGx is available
  if (!requireNamespace("PharmacoGx", quietly = TRUE)) {
    stop("PharmacoGx package required.")
  }
  
  # check if psetName is a valid name for an existing PSet
  availablePSets <- PharmacoGx::availablePSets()
  
  
  # handle error for wrong PSet name
  valid_names <- availablePSets[["PSet Name"]]
  if (!(psetName %in% valid_names)) {
    message(psetName, "' is not a valid PSet name.\n")
    message("Here are the available PharmacoSets from PharmacoGx:\n")
    print(valid_names)
    stop("\nPlease choose one of the above dataset names.")
  }
  
  
  # download requested PSet
  message("Downloading PharmacoSet: ", psetName, " ...")
  pset <- PharmacoGx::downloadPSet(psetName)
  
  # validate download was successful
  if (!inherits(pset, "PharmacoSet")) {
    stop("The downloaded object is not a valid PharmacoSet.")
  }
  
  message("Successfully loaded PharmacoSet: ", psetName)
  
  
  # Do any other preparation of the PSet for use in this package
  # TODO --> not sure yet... add here as other functions are created
  
  return(pset)
  
}