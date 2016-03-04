# Not exported
get_nuc_freqs <- function(seqs){
  #get nucleotide frequencies
  nucFreqs <- colSums(Biostrings::letterFrequency(seqs, c("A","C","G","T")))
  nucFreqs <- nucFreqs/sum(nucFreqs)
  return(nucFreqs)
}

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


#' get_motif_indices
#' 
#' Function to get indices of peaks that contain motifs
#' @param motifs \code{\link[TFBSTools]{PFMatrixList}} or \code{\link[TFBSTools]{PWMatrixList}}
#' @param peaks \code{\link[GenomicRanges]{GenomicRanges}} or \code{\link{fragmentCounts}}
#' @param genome \code{\link[BSgenome]{BSgenome}} object
#' @param p.cutoff default is 0.00005
#' @return A list with a vector of peak indices for each motif.
#'@export
get_motif_indices <- function(motifs, 
                               peaks, 
                               genome =  BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19, 
                               p.cutoff = 0.00005){
  if (inherits(peaks,"fragmentCounts")){
    peaks <- peaks@peaks
  }
  if (!inherits(peaks,"GenomicRanges")){
    stop("peaks input must be fragmentCounts or GenomicRanges.")
  }
  # If PFMatrixList or PFMatrix, convert to PWM
  if (inherits(motifs,"PFMatrixList")){
    motifs <- do.call(TFBSTools::PWMatrixList,lapply(motifs, TFBSTools::toPWM))
  } else if (inherits(motifs, "PFMatrix")){
    motifs <- TFBSTools::toPWM(motifs)
  }
  # Check class
  if (!(inherits(motifs, "PWMatrixList") | inherits(motifs, "PWMatrix"))){
    stop("motifs must be PFMatrixList, PWMatrixList, PFMatrix, or PWMatrix object")
  }
  
  seqs <- Biostrings::getSeq(genome, peaks)
  
  #get nucleotide frequencies
  nucFreqs <-  get_nuc_freqs(seqs)
  
  seqs <- as.character(seqs)
  
  motif_mats <- lapply(motifs, function(x) TFBSTools::as.matrix(x))
  
  if (is.installed("BiocParallel")){

    motif_ix = BiocParallel::bplapply(motif_mats,
                                      motif_match,
                                      seqs,
                                      nucFreqs,
                                      p.cutoff)
  } else{
    motif_ix <- lapply(motif_mats,
                       motif_match,
                       seqs,
                       nucFreqs,
                       p.cutoff)
  }
  
  return(motif_ix)
}
