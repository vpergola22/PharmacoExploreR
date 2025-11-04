# Load a PharmacoSet object from PharmacoGx
#
# @description TODO
#
#
# @param psetName the name of a PSet from PharmacoDB indicating the PSet to be 
#   loaded. This should match a value in the "PSet Name" column of the data
#   frame produced by availablePSets().
#
# @return TODO
# 
#
# @examples TODO
#
#
# @references
# Smirnov, Petr, et al. PharmacoDB: An Integrative Database for Mining in Vitro 
# Anticancer Drug Screening Studies. Vol. 46, no. D1, 9 Oct. 2017, pp. 
# D994–D1002, https://doi.org/10.1093/nar/gkx911.
# ---. “PharmacoGx: An R Package for Analysis of Large Pharmacogenomic 
# Datasets.” Bioinformatics, vol. 32, no. 8, 9 Dec. 2015, pp. 1244–1246, 
# https://doi.org/10.1093/bioinformatics/btv723.
# TODO check if this is the correct citation format
#
# @export
# @import PharmacoGx
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