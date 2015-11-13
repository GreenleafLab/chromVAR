


# kmerMapping <- setClass("kmerMapping",
#                             slots = c(mapping = 'Matrix::dgCMatrix',
#                                       colMapping = 'vector',
#                                       k = 'numeric',
#                                       m = 'numeric',
#                                       l = 'numeric'))
#
# setMethod("show",
#           signature="kmerMapping",
#           definition = function(object){
#             cat("An object of class ", class(object), "\n", sep = "")
#             cat("Mapping of kmers with k = ", object@k, " to kmers of the same length with up to m = ",
#                 object@m, " possible consecutive mismatches. \n", sep = "")
#             invisible(NULL)
#           })
#
#
# setMethod("[", signature = signature(x = "kmerMapping", i = "character", j = "missing"),
#           definition = function(x, i, j) {
#             x@mapping[i, colnames(x@mapping)]
#             })
#
# setMethod("[", signature = signature(x = "kmerMapping", i = "character", j = "character"),
#           definition = function(x, i ,j ) {
#             stopifnot(is.character(i), is.character(j))
#             x@mapping[i, x@colMapping[j]]
#           })
#
# setMethod("[", signature = signature(x = "kmerMapping", i = "DNAStringSet", j = "character"),
#           definition = function(x, i, j ) {
#             x@mapping[as.character(i), x@colMapping[j]]
#           })
#
# setMethod("[", signature = signature(x = "kmerMapping", i = "DNAStringSet", j = "missing"),
#           definition = function(x, i ,j ) {
#             x@mapping[as.character(i), colnames(x@mapping)]
#           })
#
# setMethod("[", signature = signature(x = "kmerMapping", i = "DNAStringSet", j = "DNAStringSet"),
#           definition = function(x, i ,j) {
#             x@mapping[as.character(i), x@colMapping[as.character(j)]]
#           })
#
# setMethod("[", signature = signature(x = "kmerMapping", i = "character", j = "DNAStringSet"),
#           definition = function(x, i, j) {
#             x@mapping[i, x@colMapping[as.character(j)]]
#           })
#
# setMethod("[", signature = signature(x = "kmerMapping", i = "ANY"),
#           definition = function(x, i) {
#            stop(paste("Must supply list of kmers as first argument.",
#                       " kmers list must be either character vector or DNAStringSet", sep="\n", collapse=""), call. = F)
#           })
#
#












