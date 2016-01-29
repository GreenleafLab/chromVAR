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
#' @param BPPARAM parameter for \code{\link[BiocParallel]{BiocParallel}}
#' @param p.cutoff default is 0.00005
#' @return A list with a vector of peak indices for each motif.
#' @export
get_motif_indices <- function(motifs, 
                              peaks, 
                              genome =  BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19, 
                              BPPARAM = BiocParallel::bpparam(), 
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
  seq_widths <- Biostrings::width(seqs)
  n_seq <- length(seqs)
  ranges <- IRanges::IRanges(start = cumsum(c(1,seq_widths)[1:n_seq]), 
                               width = seq_widths)
  tmpseq <- paste(seqs, collapse="")
  
  #get nucleotide frequencies
  nucFreqs <-  get_nuc_freqs(seqs)
  
  get_motif_matches <- function(motif, tmpseq, ranges, nucFreqs){    
    #find motifs on forward and reverse strands
    motif = as.matrix(motif)
    min.score = TFMPvalue::TFMpv2sc(motif, p.cutoff, bg = nucFreqs, type="PWM")
    forward <- Biostrings::matchPWM(motif, tmpseq,  
                                    min.score = min.score)   
    forward_matches <- IRanges::findOverlaps(forward,ranges,
                                             type="within",
                                             select="all")
    reverse <- Biostrings::matchPWM(Biostrings::reverseComplement(motif),
                                    tmpseq, 
                                    min.score = min.score)    
    reverse_matches <- IRanges::findOverlaps(reverse,ranges,
                                             type="within",
                                             select="all")

    out <- unique(c(IRanges::subjectHits(forward_matches),
                    IRanges::subjectHits(reverse_matches)))
    return(out)
  }  
    
  #find motif matches
  if (is.installed("BiocParallel")){
    motif_ix <- BiocParallel::bplapply(motifs,
                                     get_motif_matches,
                                     tmpseq,
                                     ranges,
                                     nucFreqs,
                                     BPPARAM = BPPARAM)    
  } else{
    motif_ix <- lapply(motifs,
                       get_motif_matches,
                       tmpseq,
                       ranges,
                       nucFreqs,
                       BPPARAM = BPPARAM)  
  }
  
  return(motif_ix)
}






