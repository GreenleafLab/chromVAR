#' getJasparMotifs
#'
#' Function to get motifs from JASPAR database
#' @param species Which species?  use eithe jaspar code or latin name.
#' default is 'Homo sapiens'
#' @param collection Which collection to use?  default is 'CORE'
#' @param jaspar_db Which JASPAR database? default is `JASPAR2016::JASPAR2016`
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
#' # use the JASPAR2020 database
#' motifs <- getJasparMotifs(jaspar_db = JASPAR2020::JASPAR2020)
#'
getJasparMotifs <- function(species = "Homo sapiens",
                            collection = "CORE",
                            jaspar_db = JASPAR2016::JASPAR2016,
                            ...) {
  opts <- list()
  opts["species"] <- species
  opts["collection"] <- collection
  opts <- c(opts, list(...))
  out <- TFBSTools::getMatrixSet(jaspar_db, opts)
  if (!isTRUE(all.equal(TFBSTools::name(out), names(out))))
    names(out) <- paste(names(out), TFBSTools::name(out), sep = "_")
  return(out)
}
