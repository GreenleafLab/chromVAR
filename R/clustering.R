#' getSampleDistance
#'
#' Get distance between samples based on bias corrected deviations
#' @param object deviations result
#' @param threshold threshold for variability
#' @param initial_dims initial dimentions for preliminary dimensionality 
#' reduction via pca
#' @param distance_function distance function to use
#' @details This function will compute the distance between samples based on the
#'  normalized deviations.  It will first remove correlated motifs / peak sets. 
#' Then the dimensionality will be further reduced via PCA if the number of 
#' dimensions exceeds initial_dims.  Then  the supplied distance_function will 
#' be used.  
#' @return dist object for distance between samples
#' @export
#' @author Alicia Schep
#' @seealso \code{\link{getSampleCorrelation}}
#' @examples 
#' # Load very small example results from computeDeviations
#' data(mini_dev, package = "chromVAR")
#' sample_dist <- getSampleDistance(mini_dev, threshold = 0.8)  
#' # setting very low variabilitiy threshold because this is mini data set
#' # threshold should generally be above 1
#' # Use plotVariability to get a sense of an appropriate threshold
#' # As this is mini data set, results  not meaningful!
getSampleDistance <- function(object, 
                                threshold = 1.5, 
                                initial_dims = 50, 
                                distance_function = dist) {
  stopifnot(is(object, "chromVARDeviations") || 
              canCoerce(object, "chromVARDeviations"))
  stopifnot(initial_dims >= 1)
  vars <- row_sds(deviationScores(object), na_rm = TRUE)
  ix <- which(vars >= threshold)
  ix2 <- ix[remove_correlated_helper(deviations(object)[ix, , drop = FALSE], 
                                     vars[ix])]
  if (initial_dims < length(ix2)) {
    pc_res <- prcomp(t(deviations(object)[ix2, ]))
    mat <- pc_res$x[, seq_len(initial_dims)]
  } else {
    mat <- t(deviations(object)[ix2, , drop = FALSE])
  }
  d <- distance_function(mat)
  return(d)
}


#' getSampleCorrelation
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
#' @seealso \code{\link{getSampleDistance}}
#' @examples 
#' # Load very small example results from computeDeviations
#' data(mini_dev, package = "chromVAR")
#' sample_cor <- getSampleCorrelation(mini_dev, threshold = 0.8)  
#' # setting very low variabilitiy threshold because this is mini data set
#' # threshold should generally be above 1
#' # Use plotVariability to get a sense of an appropriate threshold
#' # As this is mini data set, results probably not meaningful!
getSampleCorrelation <- function(object, threshold = 1.5) {
  stopifnot(is(object, "chromVARDeviations") || 
              canCoerce(object, "chromVARDeviations"))
  vars <- row_sds(deviationScores(object), na_rm = TRUE)
  ix <- which(vars >= threshold)
  ix2 <- ix[remove_correlated_helper(deviations(object)[ix, , drop = FALSE], 
                                     vars[ix])]
  cormat <- cor(deviations(object)[ix2, , drop = FALSE], 
                use = "pairwise.complete.obs")
  return(cormat)
}
