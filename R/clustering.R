
get_sample_distance <- function(object, threshold = 1.5, initial_dims = 50, distance_function = dist){
  vars <- chromVAR:::row_sds(assays(object)$z)
  ix <- which(vars >= threshold)
  ix2 <- ix[chromVAR:::remove_correlated_helper(assays(object)$deviations[ix,], vars[ix])]
  if (initial_dims < length(ix2)){
    pc_res <- prcomp(t(assays(object)$deviations[ix2,]))
    mat <- pc_res$x[,1:initial_dims]
  } else{
    mat <- t(assays(object)$deviations[ix2,])
  }
  d <- distance_function(mat)
  return(d)
}


get_sample_correlation <- function(object, threshold = 1.5){
  vars <- chromVAR:::row_sds(assays(object)$z)
  ix <- which(vars >= threshold)
  ix2 <- ix[chromVAR:::remove_correlated_helper(assays(object)$deviations[ix,], vars[ix])]
  cormat <- cor(assays(object)$deviations[ix2,], use = "pairwise.complete.obs")
  return(cormat)
}
