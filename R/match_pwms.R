
#' match_pwms
#'
#' Find pwm matches
#' @param pwms either \code{\link[TFBSTools]{PFMatrix}}, \code{\link[TFBSTools]{PFMatrixList}},
#' \code{\link[TFBSTools]{PWMatrix}}, \code{\link[TFBSTools]{PWMatrixList}}
#' @param subject either \code{\link[GenomicRanges]{GenomicRanges}}, \code{\link[Biostrings]{DNAStringSet}},
#' \code{\link[Biostrings]{DNAString}}, or character vector
#' @param genome BSgenome object, only used if subect is \code{\link[GenomicRanges]{GenomicRanges}}
#' @param bg background nucleotide frequencies. if not provided, computed from subject
#' @param out what to return? see details
#' @param p.cutoff p-value cutoff for returning motifs
#' @param w parameter controlling size of window for filtration; default is 7
#' @details Can either return a sparse matrix with values set to 1 for a match (if return == "match"), a sparse matrix
#' with values set to the max motif score in each sequence (but zero for sequences with no score above minimum p value threshold),
#' or \code{\link[GenomicRanges]{GenomicRanges}} if positions
#' @import Biostrings
#' @import TFBSTools
#' @import Matrix
#'@export
setGeneric("match_pwms", function(pwms, subject,...) standardGeneric("match_pwms"))

#' @describeIn match_pwms
#' @export
setMethod("match_pwms", signature(pwms = "PWMatrixList", subject = "DNAStringSet"),
          function(pwms, subject,  bg = NULL, out = c("match","scores","positions"),p.cutoff = 0.00005, w =7){
            out = match.arg(out)

            if (is.null(bg)){
              bg <- get_nuc_freqs(subject)
            } else{
              stopifnot(length(bg)==4 && is.numeric(bg))
            }
            seqs <- as.character(subject)
            motif_mats <- convert_pwms(pwms, bg)

            if (out == "match"){
              out <- get_motif_ix(motif_mats,seqs,bg,p.cutoff,w)
              colnames(out) <- names(pwms)
            } else if (out == "scores"){
              out <- get_max_motif_score(motif_mats,seqs,bg,p.cutoff,w)
              colnames(out) <- names(pwms)
            } else{
              stop("positions not yet implemented")
            }
            return(out)
            })

#' @describeIn match_pwms
#' @export
setMethod("match_pwms", signature(pwms = "PWMatrixList", subject = "character"),
          function(pwms, subject,  bg = NULL, out = c("match","scores","positions"),p.cutoff = 0.00005, w =7){
            out = match.arg(out)

            if (is.null(bg)){
              bg <- get_nuc_freqs(subject)
            } else{
              stopifnot(length(bg)==4 && is.numeric(bg))
            }
            motif_mats <- convert_pwms(pwms, bg)

            if (out == "match"){
              out <- get_motif_ix(motif_mats,subject,bg,p.cutoff,w)
              colnames(out) <- names(pwms)
            } else if (out == "scores"){
              out <- get_max_motif_score(motif_mats,seqs,bg,p.cutoff,w)
              colnames(out) <- names(pwms)
            } else{
              stop("positions not yet implemented")
            }
            return(out)
          })

#' @describeIn match_pwms
#' @export
setMethod("match_pwms", signature(pwms = "PWMatrixList", subject = "DNAString"),
          function(pwms, subject,  bg = NULL, out = c("match","scores","positions"),p.cutoff = 0.00005, w = 7){
            out = match.arg(out)
            if (is.null(bg)){
              bg <- get_nuc_freqs(subject)
            } else{
              stopifnot(length(bg)==4 && is.numeric(bg))
            }
            motif_mats <- convert_pwms(pwms, bg)
            seqs <- as.character(subject)

            if (out == "match"){
              out <- get_motif_ix(motif_mats,seqs,bg,p.cutoff,w)
              colnames(out) <- names(pwms)
            } else if (out == "scores"){
              out <- get_max_motif_score(motif_mats,seqs,bg,p.cutoff,w)
              colnames(out) <- names(pwms)
            } else{
              stop("positions not yet implemented")
            }
            return(out)
          })

#' @describeIn match_pwms
#' @export
setMethod("match_pwms", signature(pwms = "PWMatrixList", subject = "GenomicRanges"),
          function(pwms, subject, genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19,  bg = NULL, out = c("match","scores","positions"),p.cutoff = 0.00005, w = 7){

            out = match.arg(out)
            GenomicRanges::strand(subject) <- "+"
            seqs <- getSeq(genome, subject)
            if (is.null(bg)){
              bg <- get_nuc_freqs(seqs)
            } else{
              stopifnot(length(bg)==4 && is.numeric(bg))
            }
            motif_mats <- convert_pwms(pwms, bg)
            seqs <- as.character(seqs)

            if (out == "match"){
              out <- get_motif_ix(motif_mats,seqs,bg,p.cutoff,w)
              colnames(out) <- names(pwms)
            } else if (out == "scores"){
              out <- get_max_motif_score(motif_mats,seqs,bg,p.cutoff,w)
              colnames(out) <- names(pwms)
            } else{
              tmp_out <- get_motif_positions(motif_mats,seqs,bg,p.cutoff,w)
              out <- lapply(1:length(motif_mats), function(x){
                m_ix <- which(tmp_out$motif_ix == x - 1)
                GenomicRanges::GRanges(GenomicRanges::seqnames(subject)[tmp_out$seq_ix[m_ix] + 1],
                                              IRanges::IRanges(start = GenomicRanges::start(subject[tmp_out$seq_ix[m_ix] + 1]) + tmp_out$pos[m_ix],
                                                               width = ncol(motif_mats[[x]])),
                                                  strand = tmp_out$strand[m_ix], score = tmp_out$score[m_ix])
              })
              names(out) <- names(pwms)
            }
            return(out)
          })


### PFMatrixList ---------------------------------------------------------------


#' @describeIn match_pwms
#' @export
setMethod("match_pwms", signature(pwms = "PFMatrixList", subject = "DNAStringSet"),
          function(pwms, subject, bg = NULL, out = c("match","scores","positions"), p.cutoff = 0.00005, w =7){
            out = match.arg(out)

            if (is.null(bg)){
              bg <- get_nuc_freqs(subject)
            } else{
              stopifnot(length(bg)==4 && is.numeric(bg))
            }
            motif_mats <- lapply(pwms, function(x) as.matrix(toPWM(x, bg = bg)))
            seqs <- as.character(subject)

            if (out == "match"){
              out <- get_motif_ix(motif_mats,seqs,bg,p.cutoff,w)
              colnames(out) <- names(pwms)
            } else if (out == "scores"){
              out <- get_max_motif_score(motif_mats,seqs,bg,p.cutoff,w)
              colnames(out) <- names(pwms)
            } else{
              stop("positions not yet implemented")
            }
            return(out)
          })

#' @describeIn match_pwms
#' @export
setMethod("match_pwms", signature(pwms = "PFMatrixList", subject = "character"),
          function(pwms, subject, bg = NULL, out = c("match","scores","positions"), p.cutoff = 0.00005, w =7){
            out = match.arg(out)

            if (is.null(bg)){
              bg <- get_nuc_freqs(subject)
            } else{
              stopifnot(length(bg)==4 && is.numeric(bg))
            }
            motif_mats <- lapply(pwms, function(x) as.matrix(toPWM(x, bg = bg)))
            if (out == "match"){
              out <- get_motif_ix(motif_mats,subject,bg,p.cutoff,w)
              colnames(out) <- names(pwms)
            } else if (out == "scores"){
              out <- get_max_motif_score(motif_mats,seqs,bg,p.cutoff,w)
              colnames(out) <- names(pwms)
            } else{
              stop("positions not yet implemented")
            }
            return(out)
          })

#' @describeIn match_pwms
#' @export
setMethod("match_pwms", signature(pwms = "PFMatrixList", subject = "DNAString"),
          function(pwms, subject,  bg = NULL, out = c("match","scores","positions"),p.cutoff = 0.00005, w = 7){
            out = match.arg(out)

            if (is.null(bg)){
              bg <- get_nuc_freqs(subject)
            } else{
              stopifnot(length(bg)==4 && is.numeric(bg))
            }
            seqs <- as.character(subject)

            motif_mats <- lapply(pwms, function(x) as.matrix(toPWM(x, bg = bg)))

            if (out == "match"){
              out <- get_motif_ix(motif_mats,seqs,bg,p.cutoff,w)
              colnames(out) <- names(pwms)
            } else if (out == "scores"){
              out <- get_max_motif_score(motif_mats,seqs,bg,p.cutoff,w)
              colnames(out) <- names(pwms)
            } else{
              stop("positions not yet implemented")
            }
            return(out)
          })

#' @describeIn match_pwms
#' @export
setMethod("match_pwms", signature(pwms = "PFMatrixList", subject = "GenomicRanges"),
          function(pwms, subject, genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19, bg = NULL, out = c("match","scores","positions"), p.cutoff = 0.00005, w = 7){

            out = match.arg(out)
            GenomicRanges::strand(subject) = "+"
            seqs <- getSeq(genome, subject)
            if (is.null(bg)){
              bg <- get_nuc_freqs(seqs)
            } else{
              stopifnot(length(bg)==4 && is.numeric(bg))
            }
            seqs <- as.character(seqs)

            motif_mats <- lapply(pwms, function(x) as.matrix(toPWM(x, bg = bg)))
            if (out == "match"){
              out <- get_motif_ix(motif_mats,seqs,bg,p.cutoff,w)
              colnames(out) <- names(pwms)
            } else if (out == "scores"){
              out <- get_max_motif_score(motif_mats,seqs,bg,p.cutoff,w)
              colnames(out) <- names(pwms)
            } else{
              tmp_out <- get_motif_positions(motif_mats,seqs,bg,p.cutoff,w)
              out <- lapply(1:length(motif_mats), function(x){
                m_ix <- which(tmp_out$motif_ix == x - 1)
                GenomicRanges::GRanges(GenomicRanges::seqnames(subject)[tmp_out$seq_ix[m_ix] + 1],
                                                IRanges::IRanges(start = GenomicRanges::start(subject[tmp_out$seq_ix[m_ix] + 1]) + tmp_out$pos[m_ix],
                                                                 width = ncol(motif_mats[[x]])),
                                                strand = tmp_out$strand[m_ix], score = tmp_out$score[m_ix])

              })
              names(out) <- names(pwms)
            }
            return(out)
          })

# Single PWM input -------------------------------------------------------------

#' @describeIn match_pwms
#' @export
setMethod("match_pwms", signature(pwms = "PWMatrix", subject = "DNAStringSet"),
          function(pwms, subject, bg = NULL, out = c("match","scores","positions"), p.cutoff = 0.00005, w = 7){
            out = match.arg(out)
            if (is.null(bg)){
              bg <- get_nuc_freqs(subject)
            } else{
              stopifnot(length(bg)==4 && is.numeric(bg))
            }
            motif_mats <- convert_pwms(pwms, bg)
            seqs <- as.character(subject)

            if (out == "match"){
              out <- get_motif_ix(motif_mats,seqs,bg,p.cutoff,w)
            } else if (out == "scores"){
              out <- get_max_motif_score(motif_mats,seqs,bg,p.cutoff,w)
            } else{
              stop("positions not yet implemented")
            }
            return(out)
          })

#' @describeIn match_pwms
#' @export
setMethod("match_pwms", signature(pwms = "PWMatrix", subject = "character"),
          function(pwms, subject, bg = NULL, out = c("match","scores","positions"), p.cutoff = 0.00005, w =7){
            out = match.arg(out)

            if (is.null(bg)){
              bg <- get_nuc_freqs(subject)
            } else{
              stopifnot(length(bg)==4 && is.numeric(bg))
            }

            motif_mats <- convert_pwms(pwms, bg)

            if (out == "match"){
              out <- get_motif_ix(motif_mats,subject,bg,p.cutoff,w)
            } else if (out == "scores"){
              out <- get_max_motif_score(motif_mats,subject,bg,p.cutoff,w)
            } else{
              stop("positions not yet implemented")
            }
            return(out)
          })

#' @describeIn match_pwms
#' @export
setMethod("match_pwms", signature(pwms = "PWMatrix", subject = "DNAString"),
          function(pwms, subject, bg = NULL, out = c("match","scores","positions"), p.cutoff = 0.00005, w = 7){
            out = match.arg(out)

            if (is.null(bg)){
              bg <- get_nuc_freqs(subject)
            } else{
              stopifnot(length(bg)==4 && is.numeric(bg))
            }
            seqs <- as.character(subject)
            motif_mats <- convert_pwms(pwms, bg)

            if (out == "match"){
              out <- get_motif_ix(motif_mats,seqs,bg,p.cutoff,w)
            } else if (out == "scores"){
              out <- get_max_motif_score(motif_mats,seqs,bg,p.cutoff,w)
            } else{
              stop("positions not yet implemented")
            }
            return(out)
          })

#' @describeIn match_pwms
#' @export
setMethod("match_pwms", signature(pwms = "PWMatrix", subject = "GenomicRanges"),
          function(pwms, subject, genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19, bg = NULL, out = c("match","scores","positions"), p.cutoff = 0.00005, w = 7){

            out = match.arg(out)
            GenomicRanges::strand(subject) <- "+"
            seqs <- getSeq(genome, subject)
            if (is.null(bg)){
              bg <- get_nuc_freqs(seqs)
            } else{
              stopifnot(length(bg)==4 && is.numeric(bg))
            }
            seqs <- as.character(seqs)
            motif_mats <- convert_pwms(pwms, bg)

            if (out == "match"){
              out <- get_motif_ix(motif_mats,seqs,bg,p.cutoff,w)
            } else if (out == "scores"){
              out <- get_max_motif_score(motif_mats,seqs,bg,p.cutoff,w)
            } else{
              tmp_out <- get_motif_positions(motif_mats,seqs,bg,p.cutoff,w)
              m_ix <- which(tmp_out$motif_ix == x - 1)
              out <- GenomicRanges::GRanges(GenomicRanges::seqnames(subject)[tmp_out$seq_ix[m_ix] + 1],
                                              IRanges::IRanges(start = GenomicRanges::start(subject[tmp_out$seq_ix[m_ix] + 1]) + tmp_out$pos[m_ix],
                                                               width = ncol(motif_mats[[x]])),
                                              strand = tmp_out$strand[m_ix], score = tmp_out$score[m_ix])

            }
            return(out)
          })

# Single PFM -------------------------------------------------------------------

#' @describeIn match_pwms
#' @export
setMethod("match_pwms", signature(pwms = "PFMatrix", subject = "DNAStringSet"),
          function(pwms, subject, bg = NULL, out = c("match","scores","positions"), p.cutoff = 0.00005, w = 7){
            out = match.arg(out)

            if (is.null(bg)){
              bg <- get_nuc_freqs(subject)
            } else{
              stopifnot(length(bg)==4 && is.numeric(bg))
            }
            motif_mats <- list(as.matrix(toPWM(pwm, bg = bg)))
            seqs <- as.character(subject)

            if (out == "match"){
              out <- get_motif_ix(motif_mats,seqs,bg,p.cutoff,w)
            } else if (out == "scores"){
              out <- get_max_motif_score(motif_mats,seqs,bg,p.cutoff,w)
            } else{
              stop("positions not yet implemented")
            }
            return(out)
          })

#' @describeIn match_pwms
#' @export
setMethod("match_pwms", signature(pwms = "PFMatrix", subject = "character"),
          function(pwms, subject, bg = NULL, out = c("match","scores","positions"), p.cutoff = 0.00005, w =7){
            out = match.arg(out)

            if (is.null(bg)){
              bg <- get_nuc_freqs(subject)
            } else{
              stopifnot(length(bg)==4 && is.numeric(bg))
            }
            motif_mats <- list(as.matrix(toPWM(pwm, bg = bg)))
            if (out == "match"){
              out <- get_motif_ix(motif_mats,subject,bg,p.cutoff,w)
            } else if (out == "scores"){
              out <- get_max_motif_score(motif_mats,subject,bg,p.cutoff,w)
            } else{
              stop("positions not yet implemented")
            }
            return(out)
          })

#' @describeIn match_pwms
#' @export
setMethod("match_pwms", signature(pwms = "PFMatrix", subject = "DNAString"),
          function(pwms, subject, bg = NULL, out = c("match","scores","positions"), p.cutoff = 0.00005, w = 7){
            out = match.arg(out)

            if (is.null(bg)){
              bg <- get_nuc_freqs(subject)
            } else{
              stopifnot(length(bg)==4 && is.numeric(bg))
            }
            seqs <- as.character(subject)

            motif_mats <- list(as.matrix(toPWM(pwm, bg = bg)))

            if (out == "match"){
              out <- get_motif_ix(motif_mats,seqs,bg,p.cutoff,w)
            } else if (out == "scores"){
              out <- get_max_motif_score(motif_mats,seqs,bg,p.cutoff,w)
            } else{
              stop("positions not yet implemented")
            }
            return(out)
          })

#' @describeIn match_pwms
#' @export
setMethod("match_pwms", signature(pwms = "PFMatrix", subject = "GenomicRanges"),
          function(pwms, subject, genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19, bg = NULL, out = c("match","scores","positions"), p.cutoff = 0.00005, w = 7){

            out = match.arg(out)
            GenomicRanges::strand(subject) = "+"
            seqs <- getSeq(genome, subject)
            if (is.null(bg)){
              bg <- get_nuc_freqs(seqs)
            } else{
              stopifnot(length(bg)==4 && is.numeric(bg))
            }
            seqs <- as.character(seqs)

            motif_mats <- list(as.matrix(toPWM(pwm, bg = bg)))
            if (out == "match"){
              out <- get_motif_ix(motif_mats,seqs,bg,p.cutoff,w)
            } else if (out == "scores"){
              out <- get_max_motif_score(motif_mats,seqs,bg,p.cutoff,w)
            } else{
              tmp_out <- get_motif_positions(motif_mats,seqs,bg,p.cutoff,w)
              m_ix <- which(tmp_out$motif_ix == x - 1)
              out <- GenomicRanges::GRanges(GenomicRanges::seqnames(subject)[tmp_out$seq_ix[m_ix] + 1],
                                              IRanges::IRanges(start = GenomicRanges::start(subject[tmp_out$seq_ix[m_ix] + 1]) + tmp_out$pos[m_ix],
                                                               width = ncol(motif_mats[[x]])),
                                              strand = tmp_out$strand[m_ix], score = tmp_out$score[m_ix])
            }
            return(out)
          })


