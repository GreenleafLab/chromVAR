#' get_sample_distance
#'
#' Get distance between samples based on bias corrected deviations
#' @param object deviations result
#' @param threshold threshold for variability
#' @param initial_dims initial dimentions for preliminary dimensionality 
#' reduction via pca
#' @param distance_function distance function to use
#' @details This function will compute the distance between samples based on the 
#' normalized deviations.  It will first remove correlated motifs / peak sets. 
#' Then the dimensionality will be further reduced via PCA if the number of 
#' dimensions exceeds initial_dims.  Then  the supplied distance_function will 
#' be used.  
#' @return dist object for distance between samples
#' @export
#' @author Alicia Schep
#' @seealso \code{\link{get_sample_correlation}}
#' @examples 
#' # Load very small example counts (already filtered)
#' data(mini_counts, package = "chromVAR")
#' motifs <- get_jaspar_motifs()[c(1,2,4,298)] # only use a few for demo 
#' library(motifmatchr)
#' motif_ix <- match_motifs(motifs, mini_counts)
#'
#' # computing deviations
#' dev <- compute_deviations(object = mini_counts, 
#'                          annotations = motif_ix)
#' sample_dist <- get_sample_distance(dev, threshold = 1)  
#' # setting very low variabilitiy threshold because this is mini data set
#' # Use plot_variability to get a sense of an appropriate threshold
#' # As this is mini data set, results probably not very meaningful!
get_sample_distance <- function(object, 
                                threshold = 1.5, 
                                initial_dims = 50, 
                                distance_function = dist) {
  stopifnot(is(object, "chromVARDeviations") || 
              canCoerce(object, "chromVARDeviations"))
  stopifnot(initial_dims >= 1)
  vars <- row_sds(deviation_scores(object), na_rm = TRUE)
  ix <- which(vars >= threshold)
  ix2 <- ix[remove_correlated_helper(deviations(object)[ix, , drop = FALSE], 
                                     vars[ix])]
  if (initial_dims < length(ix2)) {
    pc_res <- prcomp(t(deviations(object)[ix2, ]))
    mat <- pc_res$x[, 1:initial_dims]
  } else {
    mat <- t(deviations(object)[ix2, , drop = FALSE])
  }
  d <- distance_function(mat)
  return(d)
}


#' get_sample_correlation
#'
#' Get correlation between samples based on bias corrected deviations
#' @param object deviations result
#' @param threshold threshold for variability
#' @details This function will compute the correlation between samples based on 
#' the normalized deviations. It will first remove correlated motifs/peak sets. 
#' Then the pearson correlation coefficient will be computed and returned.
#' @return correlation matrix between samples
#' @export
#' @author Alicia Schep
#' @seealso \code{\link{get_sample_distance}}
#' @examples 
#' # Load very small example counts (already filtered)
#' data(mini_counts, package = "chromVAR")
#' motifs <- get_jaspar_motifs()[c(1,2,4,298)] # only use a few for demo 
#' library(motifmatchr)
#' motif_ix <- match_motifs(motifs, mini_counts)
#'
#' # computing deviations
#' dev <- compute_deviations(object = mini_counts, 
#'                          annotations = motif_ix)
#' sample_cor <- get_sample_correlation(dev, threshold = 1)  
#' # setting very low variabilitiy threshold because this is mini data set
#' # Use plot_variability to get a sense of an appropriate threshold
#' # As this is mini data set, results probably not very meaningful!
get_sample_correlation <- function(object, threshold = 1.5) {
  stopifnot(is(object, "chromVARDeviations") || 
              canCoerce(object, "chromVARDeviations"))
  vars <- row_sds(deviation_scores(object), na_rm = TRUE)
  ix <- which(vars >= threshold)
  ix2 <- ix[remove_correlated_helper(deviations(object)[ix, , drop = FALSE], 
                                     vars[ix])]
  cormat <- cor(deviations(object)[ix2, , drop = FALSE], 
                use = "pairwise.complete.obs")
  return(cormat)
}
