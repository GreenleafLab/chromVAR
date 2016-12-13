


#' get_sample_distance
#'
#' @param object deviations result
#' @param threshold threshold for variability
#' @param initial_dims initial dimentions for preliminary dimensionality reduction via pca
#' @param distance_function distance function to use
#' @details This function will compute the distance between samples based on the 
#' normalized deviations.  It will first remove correlated motifs / peak sets.  Then
#' the dimensionality will be further reduced via PCA if the number of dimensions
#' exceeds initial_dims.  Then  the supplied distance_function will be use.  
#' @return dist object for distance between samples
#' @export
#'
#' @examples
get_sample_distance <- function(object, threshold = 1.5, initial_dims = 50, distance_function = dist){
  vars <- chromVAR:::row_sds(assays(object)$z)
  ix <- which(vars >= threshold)
  ix2 <- ix[remove_correlated_helper(assays(object)$deviations[ix,], vars[ix])]
  if (initial_dims < length(ix2)){
    pc_res <- prcomp(t(assays(object)$deviations[ix2,]))
    mat <- pc_res$x[,1:initial_dims]
  } else{
    mat <- t(assays(object)$deviations[ix2,])
  }
  d <- distance_function(mat)
  return(d)
}


#' get_sample_correlation
#'
#' @param object deviations result
#' @param threshold threshold for variability
#' @details This function will compute the correlation between samples based on the 
#' normalized deviations.  It will first remove correlated motifs / peak sets. Then
#' the pearson correlation coefficient will be computed and returned.
#' @return correlation matrix between samples
#' @export
#'
#' @examples
get_sample_correlation <- function(object, threshold = 1.5){
  vars <- chromVAR:::row_sds(assays(object)$z)
  ix <- which(vars >= threshold)
  ix2 <- ix[remove_correlated_helper(assays(object)$deviations[ix,], vars[ix])]
  cormat <- cor(assays(object)$deviations[ix2,], use = "pairwise.complete.obs")
  return(cormat)
}
