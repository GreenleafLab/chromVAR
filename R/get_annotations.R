#' getAnnotations
#' 
#' @param annotations matrix, Matrix, or data.frame of fragment counts,
#' or SummarizedExperiment with counts assays, see details
#' @param ... additional arguments to pass to SummarizedExperiment
#' @return SummarizedExperiment object with 'matches' assay
#' @export
#' @author Alicia Schep
#' @examples
#' 
#' # First get example counts
#' data(mini_counts, package = "chromVAR")
#' 
#' # Get annotations from genomic ranges list
#' library(GenomicRanges)
#' library(SummarizedExperiment)
#' my_annotation_granges <- GRangesList(GRanges("chr1", 
#'                                              ranges = IRanges(start = 
#'                                              c(566763,805090), width = 8)),
#'                                      GRanges("chr1", ranges = IRanges(start =
#'                                                c(566792,895798), width = 8)))
#' anno_ix <- getAnnotations(my_annotation_granges, 
#'                            rowRanges = rowRanges(mini_counts))
setGeneric("getAnnotations", 
           function(annotations, ...) standardGeneric("getAnnotations"))

#' annotationMatches
#' 
#' @param object SummarizedExperiment with matches slot, see details
#' @details Will extract matrix from the "matches", "annotationMatches", or
#' "motif_matches" assay of a SummarizedExperiment
#' @return logical matrix of annotation matches
#' @export
#' @author Alicia Schep
#' @rdname annotationMatches
#' @name annotationMatches
#' @aliases annotationMatches,SummarizedExperiment-method 
#' annoation_matches<-,SummarizedExperiment-method
#' @examples 
#' # load annotation matrix; result from matchMotifs
#' data(mini_ix, package = "chromVAR")
#' matches <- annotationMatches(mini_ix)
setGeneric("annotationMatches", 
           function(object) standardGeneric("annotationMatches"))

#' @rdname annotationMatches
setGeneric("annotationMatches<-", 
           function(object, value) standardGeneric("annotationMatches<-"))


#' @rdname annotationMatches
setMethod("annotationMatches", 
          c(object = "SummarizedExperiment"), 
          function(object) {
            #object <- matches_check(object)
            if ("annotationMatches" %in% assayNames(object)){
              out <- assays(object)$annotationMatches
            } else if ("annotation_matches" %in% assayNames(object)){
              out <- assays(object)$annotation_matches
            } else if ("motif_matches" %in% assayNames(object)){
              out <- assays(object)$motif_matches
              warning("motif_matches assay deprecated; update motifmatchr")
            } else if ("motifMatches" %in% assayNames(object)){
              out <- assays(object)$motifMatches
            } else if ("matches" %in% assayNames(object)){
              out <- assays(object)$matches
            } else {
              stop("No appropriately named assay. See Details section in man", 
                   "page")
            }
            return(out)
          })

#' @rdname annotationMatches
#' @param value logical Matrix with annotation matches
#' @export
setReplaceMethod("annotationMatches", 
          c(object = "SummarizedExperiment"), 
          function(object, value) {
            #object <- matches_check(object)
            stopifnot(canCoerce(value,"lMatrix"))
            if (!is(value, "lMatrix")) {
              value <- as(value, "lMatrix")
              warning("Annotation object matches converted to logical")
            }
            if ("annotationMatches" %in% assayNames(object)){
              assays(object)$annotationMatches <- value
            } else if ("annotation_matches" %in% assayNames(object)){
              assays(object)$annotation_matches <- value
            } else if ("motif_matches" %in% assayNames(object)){
              assays(object)$motif_matches <- value
              warning("motif_matches assay deprecated; update motifmatchr")
            } else if ("motifMatches" %in% assayNames(object)){
              assays(object)$motifMatches <- value
            }  else if ("matches" %in% assayNames(object)){
              assays(object)$matches <- value
            } else {
              assays(object)$annotationMatches <- value
            }
            return(object)
          })

#' @param rowRanges GenomicRanges or GenomicRangesList or 
#' RangedSummarizedExperiment
#' @describeIn getAnnotations get annotation matrix from GRangesList
#' @export
setMethod(getAnnotations, 
          c(annotations = "GRangesList"), function(annotations, 
                                                   rowRanges, ...) {
            if (is(rowRanges, "RangedSummarizedExperiment")) 
              rowRanges <- rowRanges(rowRanges)
            matches <- vapply(annotations, 
                              function(x) overlapsAny(rowRanges, x),
                              rep(TRUE,length(rowRanges)))
            SummarizedExperiment(assays = 
                                   list(annotationMatches = Matrix(matches)),
                                 rowRanges = rowRanges, ...)
          })

#' @describeIn getAnnotations get annotation matrix from Matrix or matrix
#' @export
setMethod(getAnnotations, c(annotations = "MatrixOrmatrix"), 
          function(annotations,  ...) {
            SummarizedExperiment(assays = 
                                   list(annotationMatches = as(annotations,
                                                                "lMatrix")),
                                 ...)
          })

#' @describeIn getAnnotations get annotation matrix from data.frame
#' @export
setMethod(getAnnotations, c(annotations = "data.frame"), 
          function(annotations, ...) {
            getAnnotations(as.matrix(annotations), ...)
          })

#' @param npeaks number of peaks
#' @describeIn getAnnotations get annotation matrix from list
#' @export
setMethod(getAnnotations, c(annotations = "list"), 
          function(annotations, npeaks = NULL, 
                   ...) {
            add_args <- list(...)
            if (is.null(npeaks) && "rowData" %ni% names(add_args) && 
                "rowRanges" %ni% names(add_args)) {
              stop("Must provide npeaks, rowData, or rowRanges")
            } else if (is.null(npeaks) && "rowRanges" %in% names(add_args)) {
              npeaks <- length(add_args[["rowRanges"]])
            } else if (is.null(npeaks)) {
              npeaks <- nrow(add_args[["rowData"]])
            }
            SummarizedExperiment(assays = 
                                   list(annotationMatches = 
                                          convert_from_ix_list(annotations, 
                                                               npeaks)),
                                 ...)
          })


#' @param column column of bed file with annotation names, see details
#' @describeIn getAnnotations get annotations from bed files
#' @export
setMethod(getAnnotations, c(annotations = "character"), 
          function(annotations, rowRanges, 
                   column = NULL, ...) {
            if (length(annotations) == 1 && !is.null(column)) {
              if (is.installed("readr")) {
                bed <- 
                  as.data.frame(suppressMessages(
                    readr::read_tsv(file = 
                                      annotations,
                                    col_names = FALSE)[, c(1:3, column)]))
              } else {
                bed <- read.delim(file = annotations, header = FALSE, 
                                  sep = "\t", 
                                  stringsAsFactors = FALSE)[, c(1:3, column)]
              }
              if (!is.null(column)) {
                colnames(bed) <- c("chr", "start", "end", "group")
                bed[, "start"] <- bed[, "start"] + 1
                grl <- split(makeGRangesFromDataFrame(bed), bed$group)
              } else {
                colnames(bed) <- c("chr", "start", "end")
                bed[, "start"] <- bed[, "start"] + 1
                grl <- GRangesList(makeGRangesFromDataFrame(bed))
              }
            } else {
              grl <- do.call(GRangesList, lapply(annotations, 
                                                 function(filename) {
                if (is.installed("readr")) {
                  bed <- as.data.frame(suppressMessages(
                    readr::read_tsv(file = filename, 
                                    col_names = FALSE)[, c(1:3, column)]))
                } else {
                  bed <- read.delim(file = filename, 
                                    header = FALSE, 
                                    sep = "\t",
                                    stringsAsFactors = FALSE)[,  c(1:3, column)]
                }
                colnames(bed) <- c("chr", "start", "end")
                bed[, "start"] <- bed[, "start"] + 1
                makeGRangesFromDataFrame(bed)
              }))
            }
            getAnnotations(grl, rowRanges, ...)
          })

