#' getJasparMotifs
#'
#' Function to get motifs from JASPAR database
#' @param species Which species?  use eithe jaspar code or latin name. 
#' default is 'Homo sapiens'
#' @param collection Which collection to use?  default is 'CORE'
#' @param ... additional arguments to opts for 
#' \code{\link[TFBSTools]{getMatrixSet}}
#' @details Simply a wrapper function for \code{\link[TFBSTools]{getMatrixSet}}
#'  that calls JASPAR2016 database using \code{\link[JASPAR2016]{JASPAR2016}}
#' @return \code{\link[TFBSTools]{PFMatrixList}}
#' @export
#' @examples 
#' 
#' motifs <- getJasparMotifs()
#' 
#' 
getJasparMotifs <- function(species = "Homo sapiens", 
                              collection = "CORE", ...) {
  opts <- list()
  opts["species"] <- species
  opts["collection"] <- collection
  opts <- c(opts, list(...))
  out <- TFBSTools::getMatrixSet(JASPAR2016::JASPAR2016, opts)
  if (!isTRUE(all.equal(TFBSTools::name(out), names(out)))) 
    names(out) <- paste(names(out), TFBSTools::name(out), sep = "_")
  return(out)
}
