#' get_annotations
#' 
#' @param annotations matrix, Matrix, or data.frame of fragment counts,
#' or SummarizedExperiment with counts assays, see details
#' @param ... additional arguments
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
#' anno_ix <- get_annotations(my_annotation_granges, 
#'                            rowRanges = rowRanges(mini_counts))
setGeneric("get_annotations", 
           function(annotations, ...) standardGeneric("get_annotations"))

#' annotation_matches
#' 
#' @param object SummarizedExperiment with matches slot
#' @return logical matrix of annotation matches
#' @export
#' @author Alicia Schep
#' # Load very small example counts (already filtered)
#' data(mini_counts, package = "chromVAR")
#' motifs <- get_jaspar_motifs()[c(1,2,4,298)] # only use a few for demo 
#' library(motifmatchr)
#' motif_ix <- match_motifs(motifs, mini_counts)
#' matches <- annotation_matches(motif_ix)
setGeneric("annotation_matches", 
           function(object) standardGeneric("annotation_matches"))



#' @describeIn annotation_matches get matches assay from SummarizedExperiment
setMethod("annotation_matches", 
          c(object = "SummarizedExperiment"), 
          function(object) {
            object <- matches_check(object)
            assays(object)$matches
          })

#' @param rowRanges GenomicRanges or GenomicRangesList or 
#' RangedSummarizedExperiment
#' @describeIn get_annotations get annotation matrix from GRangesList
#' @export
setMethod(get_annotations, 
          c(annotations = "GRangesList"), function(annotations, 
                                                   rowRanges, ...) {
            if (is(rowRanges, "RangedSummarizedExperiment")) 
              rowRanges <- rowRanges(rowRanges)
            matches <- sapply(annotations, 
                              function(x) overlapsAny(rowRanges, x))
            SummarizedExperiment(assays = list(matches = Matrix(matches)),
                                 rowRanges = rowRanges)
          })

#' @describeIn get_annotations get annotation matrix from Matrix or matrix
#' @export
setMethod(get_annotations, c(annotations = "MatrixOrmatrix"), 
          function(annotations,  ...) {
            SummarizedExperiment(assays = list(matches = as(annotations,
                                                            "lMatrix")),
                                 ...)
          })

#' @describeIn get_annotations get annotation matrix from data.frame
#' @export
setMethod(get_annotations, c(annotations = "data.frame"), 
          function(annotations, ...) {
            get_annotations(as.matrix(annotations), ...)
          })

#' @param npeaks number of peaks
#' @describeIn get_annotations get annotation matrix from list
#' @export
setMethod(get_annotations, c(annotations = "list"), 
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
            SummarizedExperiment(assays = list(matches = 
                                                 convert_from_ix_list(annotations, 
                                                                      npeaks)),
                                 ...)
          })


#' @param column column of bed file with annotation names, see details
#' @param colData data on each annotation to store along with annotation matrix
#' @describeIn get_annotations get annotations from bed files
#' @export
setMethod(get_annotations, c(annotations = "character"), 
          function(annotations, rowRanges, 
                   column = NULL, colData = NULL) {
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
                  bed <- read.delim(file = filename, header = FALSE, sep = "\t", 
                                    stringsAsFactors = FALSE)[,  c(1:3, column)]
                }
                colnames(bed) <- c("chr", "start", "end")
                bed[, "start"] <- bed[, "start"] + 1
                makeGRangesFromDataFrame(bed)
              }))
            }
            get_annotations(grl, rowRanges, colData = colData)
          })
