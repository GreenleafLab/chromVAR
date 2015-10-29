



#' Compute variance across set of annotations.
#' 
#' @param counts_mat matrix with fragment counts, columns represent samples, rows represent peaks
#' @param motif_mat matrix with annotations as columns, peaks as rows
#' @param bias vector with gc-content or other bias measurement for peaks
#' @param cutoff include all annotations with value greater/less than cutoff (greater or less determined by cutoff.type)
#' @param cutoff.type 'lower' or 'upper' to indicate whether cutoff is an upper or lower bound
#' @return A deviationResultSet object
#' 
#' 
#' 

calculate_deviations <- function(counts_mat, motif_mat, bias, cutoff = 0, cutoff.type = 'lower',
                                 niterations = 100, window_size = 2500, BPPARAM=bpparam()){
  
  motif_mat[motif_mat > cutoff] = 1
  
  niterations = niterations + 1
  
  n_tf = ncol(motif_mat)

  reads_per_cell = colSums(counts_mat)
  total_reads = sum(counts_mat)

  bg = setBackgroundParameters(counts_mat = counts_mat, bias =  bias, window = window_size)
  
  get_dev_tf <- function(tf_index){
    
    tf_vec = motif_mat[,tf_index]
    tf_count = sum(tf_vec)
    
    observed = tf_vec %*% counts_mat#colSums(counts_mat[tf_vec==1,])
    expected =  sum(observed)*reads_per_cell/total_reads
    
    deviation = (observed - expected)
      
    sampled_counts = sampleBackgroundPeaks(object = bg, annotation_vector = tf_vec, counts_mat = counts_mat, 
                                           reads = sum(observed), niterations = niterations)
    expected_sampled_counts =  outer(reads_per_cell/total_reads, colSums(sampled_counts))
    
    
    result = deviationResult(sampled_deviations = as.matrix(sampled_counts[,1:(niterations-1)] - expected_sampled_counts[,1:(niterations-1)]), 
                             extra_deviations =  as.vector(sampled_counts[,niterations] - expected_sampled_counts[,niterations]), 
                             observed_deviations = as.vector(deviation),
                             tf = colnames(motif_mat)[tf_index])
    
    result = compute_z_score(result)
    result = compute_variability(result)
    
    return(result)
  }
  
  
  results = deviationResultSet(results = bplapply(1:n_tf, get_dev_tf, BPPARAM = BPPARAM), tfs = colnames(motif_mat))
  results = adjust_p_values(results)
  
  return(results)
  
}


calculate_deviations2 <- function(counts_mat, motif_mat, bg_peaks, cutoff = 0, cutoff.type = 'lower',
                                 niterations = 50, window_size = 2500, BPPARAM=bpparam()){
  
  motif_mat[motif_mat > cutoff] =1
  
  niterations = niterations + 1
  
  n_tf = ncol(motif_mat)
  n_peaks = nrow(counts_mat)
  
  reads_per_cell = colSums(counts_mat)
  total_reads = sum(counts_mat)
  
  #bg = setBackgroundParameters(counts_mat = counts_mat, bias =  bias, window = window_size)
  
  get_dev_tf <- function(tf_index){
    
    tf_vec = motif_mat[,tf_index]
    tf_count = sum(tf_vec)
    
    observed = tf_vec %*% counts_mat#colSums(counts_mat[tf_vec==1,])
    expected =  sum(observed)*reads_per_cell/total_reads
    
    deviation = (observed - expected)
    
    tf_indices = do.call(c, sapply(1:length(tf_vec), function(x) rep(x, tf_vec[x])))
    sample_mat = sparseMatrix(j = as.vector(bg_peaks[tf_indices,1:niterations]), i = rep(1:niterations, each = tf_count), x=1, dims = c(niterations, n_peaks))
    
    sampled_counts =  t(sample_mat %*% counts_mat)
    expected_sampled_counts =  outer(reads_per_cell/total_reads, colSums(sampled_counts))
    
    result = deviationResult(sampled_deviations = as.matrix(sampled_counts[,1:(niterations-1)] - expected_sampled_counts[,1:(niterations-1)]), 
                             extra_deviations =  as.vector(sampled_counts[,niterations] - expected_sampled_counts[,niterations]), 
                             observed_deviations = as.vector(deviation),
                             tf = colnames(motif_mat)[tf_index])
    
    result = compute_z_score(result)
    result = compute_variability(result)
    
    return(result)
  }
  
  
  results = deviationResultSet(results = bplapply(1:n_tf, get_dev_tf, BPPARAM = BPPARAM), tfs = colnames(motif_mat))
  results = adjust_p_values(results)
  
  return(results)
}



calculate_deviations3 <- function(counts_mat, motif_mat, bias, cutoff = 0, cutoff.type = 'lower',
                                 niterations = 100, window_size = 2500, BPPARAM=bpparam()){
  
  motif_mat[motif_mat > cutoff] =1
  
  niterations = niterations + 1
  
  n_tf = ncol(motif_mat)
  n_peaks = nrow(counts_mat)
  n_cell = ncol(counts_mat)
  
  counts = rowSums(counts_mat)
  reads_per_cell = colSums(counts_mat)
  total_reads = sum(counts_mat)
  
  #Get sorted counts
  
  counts_sort = sort(counts, index.return=T)$ix 
  counts_sort_rev = sort(counts_sort, index.return=T)$ix
  
  #Get sorted bias
  bias_sort = sort(bias, index.return=T)$ix
  bias_sort_rev = sort(bias_sort, index.return=T)$ix
  
  
  get_dev_tf <- function(tf_index){
    
    tf_vec = motif_mat[,tf_index]
    tf_count = sum(tf_vec)
    
    sampling_vec = mean_smooth(tf_vec[counts_sort],window_size)
    sampling_vec = sampling_vec[counts_sort_rev] / sum(sampling_vec)
    
    sampled_peaks = sparseMatrix(i= rep(1:niterations, each=tf_count),
                                 j = sample(1:n_peaks, size = tf_count * niterations, prob = sampling_vec, replace = TRUE),
                                 dims = c(niterations, n_peaks), 
                                 x = 1)
    
    observed = tf_vec %*% counts_mat#colSums(counts_mat[tf_vec==1,])
    expected =  sum(observed)*reads_per_cell/total_reads
    
    deviation = (observed - expected)
    
    bias_vec = mean_smooth(tf_vec[bias_sort],window_size)
    bias_vec = bias_vec[bias_sort_rev] * tf_count / sum(bias_vec)
    
    bias_peak_vec = mean_smooth(bias_vec[counts_sort],window_size)
    bias_peak_vec = bias_peak_vec[counts_sort_rev] * tf_count / sum(bias_peak_vec)
    
    bias_term = bias_vec %*% counts_mat
    
    corr_term = bias_peak_vec %*% counts_mat
    
    correction = (bias_term - corr_term) * sum(observed)/sum(bias_term)
    
    sampled_counts = t(sampled_peaks %*% counts_mat) + matrix(correction,byrow=F,nrow=n_cell,ncol=niterations)
    expected_sampled_counts =  outer(reads_per_cell/total_reads, colSums(sampled_counts))
    
    
    result = deviationResult(sampled_deviations = as.matrix(sampled_counts[,1:(niterations-1)] - expected_sampled_counts[,1:(niterations-1)]), 
                             extra_deviations =  as.vector(sampled_counts[,niterations] - expected_sampled_counts[,niterations]), 
                             observed_deviations = as.vector(deviation),
                             tf = colnames(motif_mat)[tf_index])
    
    result = compute_z_score(result)
    result = compute_variability(result)
    
    return(result)
  }
  
  
  results = deviationResultSet(results = bplapply(1:n_tf, get_dev_tf, BPPARAM = BPPARAM), tfs = colnames(motif_mat))
  results = adjust_p_values(results)
  
  return(results)
  
}


calculate_deviations4 <- function(counts_mat, motif_indices, bg_peaks, cutoff = 0, cutoff.type = 'lower',
                                  niterations = 50, window_size = 2500, BPPARAM=bpparam()){
  
  motif_mat[motif_mat > cutoff] =1
  
  niterations = niterations + 1
  
  n_tf = length(motif_indices)
  n_peaks = nrow(counts_mat)
  
  reads_per_cell = colSums(counts_mat)
  total_reads = sum(counts_mat)
  
  #bg = setBackgroundParameters(counts_mat = counts_mat, bias =  bias, window = window_size)
  
  get_dev_tf <- function(tf_index){
    
    tf_count = length(motif_indices[[tf_index]])
    tf_vec = sparseMatrix(j = motif_indices[[tf_index]], i = rep(1,tf_count), x = 1, dims = c(1, n_peaks))
    
    observed = tf_vec %*% counts_mat
    expected =  sum(observed)*reads_per_cell/total_reads
    
    deviation = (observed - expected)
    
    sample_mat = sparseMatrix(j = as.vector(bg_peaks[motif_indices[[tf_index]],1:niterations]), i = rep(1:niterations, each = tf_count), x=1, dims = c(niterations, n_peaks))
    
    sampled_counts =  t(sample_mat %*% counts_mat)
    expected_sampled_counts =  outer(reads_per_cell/total_reads, colSums(sampled_counts))
    
    result = deviationResult(sampled_deviations = as.matrix(sampled_counts[,1:(niterations-1)] - expected_sampled_counts[,1:(niterations-1)]), 
                             extra_deviations =  as.vector(sampled_counts[,niterations] - expected_sampled_counts[,niterations]), 
                             observed_deviations = as.vector(deviation),
                             tf = names(motif_indices)[tf_index])
    
    result = compute_z_score(result)
    result = compute_variability(result)
    
    return(result)
  }
  
  
  results = deviationResultSet(results = bplapply(1:n_tf, get_dev_tf, BPPARAM = BPPARAM), tfs = names(motif_indices))
  results = adjust_p_values(results)
  
  return(results)
}




