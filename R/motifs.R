#' get_motifs
#' 
#' Function to get motifs from JASPAR database
#' @param species Which species?  use eithe jaspar code or latin name. default is "Homo sapiens"
#' @param collection Which collection to use?  default is "CORE"
#' @param ... additional arguments to opts for \code{\link[TFBSTools]{getMatrixSet}}
#' @details Simply a wrapper function for \code{\link[TFBSTools]{getMatrixSet}} that calls
#' JASPAR2014 database using \code{\link[JASPAR2014]{JASPAR2014}}
#' @return \code{\link[TFBSTools]{PFMatrixList}}
#' @export
get_motifs <- function(species = "Homo sapiens", collection = "CORE", ...){
  opts = list()
  opts['species'] = species
  opts['collection'] = collection
  opts = c(opts, list(...))
  TFBSTools::getMatrixSet(JASPAR2014::JASPAR2014, opts)
}

