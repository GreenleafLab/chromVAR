# Not exported

## get_nuc_freqs function to find background sequence --------------------------

setGeneric("get_nuc_freqs", function(subject,...) standardGeneric("get_nuc_freqs"))

setMethod("get_nuc_freqs", signature(subject = "DNAStringSet"),
          function(subject) {
            nucFreqs <- colSums(Biostrings::letterFrequency(subject, c("A","C","G","T")))
            nucFreqs <- nucFreqs/sum(nucFreqs)
            return(nucFreqs)
          })

setMethod("get_nuc_freqs", signature(subject = "DNAString"),
          function(subject) {
            nucFreqs <- Biostrings::letterFrequency(subject, c("A","C","G","T"))
            nucFreqs <- nucFreqs/sum(nucFreqs)
            return(nucFreqs)
          })

setMethod("get_nuc_freqs", signature(subject = "character"),
          function(subject) {
            if (length(subject) == 1){
              return(get_nuc_freqs(DNAString(subject)))
            } else{
              return(get_nuc_freqs(DNAStringSet(subject)))
            }
          })


## convert_pwm to adjust background in pwm model ---------------------------------

setGeneric("convert_pwms", function(pwms,...) standardGeneric("convert_pwms"))

setMethod("convert_pwms", signature(pwms = "PWMatrixList"),
          function(pwms, bg_freqs) {
            lapply(pwms, function(x){
              out = TFBSTools::as.matrix(x)
              norm_mat = matrix(log(TFBSTools::bg(x)) - log(bg_freqs), nrow = 4, ncol = ncol(out), byrow = FALSE)
              return(out - norm_mat)})
          })

setMethod("convert_pwms", signature(pwms = "PWMatrix"),
          function(pwms, bg_freqs) {
            out = TFBSTools::as.matrix(pwms)
            norm_mat = matrix(log(TFBSTools::bg(pwms)) - log(bg_freqs), nrow = 4, ncol = ncol(out), byrow = FALSE)
            return(list(out - norm_mat))
})




