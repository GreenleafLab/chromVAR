
quantile_helper <- function(values, quantiles, na.rm){
  if (na.rm){
    out <- quantile(values, quantiles, na.rm = TRUE)
  } else{
    if (!all_false(is.na(values))){
      out <- rep(NA,length(quantiles))
    } else{
      out <- quantile(values, quantiles)
    }
  }
  return(out)
}

#'compute_variability
#'
#'function to compute overall variability of motif sets across samples
#'@param deviations output from \code{\link{compute_deviations}}
#'@param bootstrap_error compute bootstrap confidence interval
#'@param bootstrap_samples number of bootstrap samples to take
#'@param bootstrap_quantiles quantiles for bootstrap 
#'@param na.rm remove NAs? default is true
#'@export
compute_variability <- function(deviations, 
                                bootstrap_error = TRUE, 
                                bootstrap_samples = 1000, 
                                bootstrap_quantiles = c(0.025,0.975),
                                na.rm = TRUE){
  

  if (is.matrix(deviations)){
    sd_deviations <- row_sds(deviations, na.rm)
    
    p_sd <- pchisq((ncol(deviations)-1) * (sd_deviations**2), 
                   df = (ncol(deviations)-1), 
                   lower.tail = FALSE)
    
    p_adj <- p.adjust(p = p_sd, method = "BH")  
    
    if (bootstrap_error){
      #check some arguments
      stopifnot(bootstrap_samples > 0 && all_whole(bootstrap_samples))
      stopifnot(all_true(bootstrap_quantiles > 0))
      stopifnot(all_true(bootstrap_quantiles < 1)) 
      stopifnot(bootstrap_quantiles[2] > bootstrap_quantiles[1])
      
      boot_sd <- replicate(bootstrap_samples, row_sds_perm(deviations, na.rm))
      sd_error <- apply(boot_sd, 1, quantile_helper, quantiles = bootstrap_quantiles, na.rm = na.rm)
      
      out <- data.frame(variability = sd_deviations, 
                        bootstrap_lower_bound = sd_error[1,], 
                        bootstrap_upper_bound = sd_error[2,],
                        p_value = p_sd,
                        p_value_adj = p_adj,
                        row.names = rownames(deviations))    
    } else{
      out <- data.frame(variability = sd_deviations,
                        p_value = p_sd,
                        p_value_adj = p_adj,
                        row.names = rownames(deviations))    
    }    
  } else{
    sd_deviations <- sd(deviations, na.rm = na.rm)
    
    p_sd <- pchisq((length(deviations)-1) * (sd_deviations**2), 
                   df = (length(deviations)-1), 
                   lower.tail = FALSE)
    
    p_adj <- p_sd  
    
    if (bootstrap_error){
      #check some arguments
      stopifnot(bootstrap_samples > 0 && all_whole(bootstrap_samples))
      stopifnot(all_true(bootstrap_quantiles > 0))
      stopifnot(all_true(bootstrap_quantiles < 1)) 
      stopifnot(bootstrap_quantiles[2] > bootstrap_quantiles[1])
      
      boot_sd <- replicate(bootstrap_samples, sd(deviations[sample(1:length(deviations),length(deviations), replace = TRUE)], na.rm = na.rm))
      sd_error <- quantile_helper(boot_sd, quantiles = bootstrap_quantiles, na.rm = na.rm)
      
      out <- c(variability = sd_deviations, 
                        bootstrap_lower_bound = sd_error[1], 
                        bootstrap_upper_bound = sd_error[2],
                        p_value = p_sd,
                        p_value_adj = p_adj)    
    } else{
      out <- c(variability = sd_deviations,
                        p_value = p_sd,
                        p_value_adj = p_adj)    
    }    
  }
  return(out)
}



# Helper function --------------------------------------------------------------

compute_variability_single <- function(peak_set,
                                       counts_mat,
                                       background_peaks,
                                       expectation,
                                       counts_info,
                                       norm = TRUE){
  
  #require Matrix (for some multiprocessing options)
  suppressPackageStartupMessages(library(Matrix, quietly = TRUE, warn.conflicts = FALSE))
  
  if (length(peak_set) == 0){
    return(NA)
  }
  
  ### counts_mat should already be normed!
  tf_count <- length(peak_set)
  ### Determine if any sets have too low expected counts
  expected_totals <- sum(expectation[peak_set]) * counts_info$fragments_per_sample
  fail_filter <- which(expected_totals < 5)
  
  if (tf_count == 1){
    print('bah')
    observed <- as.vector(counts_mat[peak_set,])
    expected <- expectation[peak_set] * counts_info$fragments_per_sample
    if (norm){
      expected <- expected / sqrt(expected)
    }
    observed_deviation <- observed - expected
    sampled <- counts_mat[background_peaks[peak_set,],]
    if (norm){
      sampled_expected <-  outer(expectation[background_peaks[peak_set,]] / sqrt(expectation[background_peaks[peak_set,]]),
                                 counts_info$fragments_per_sample / sqrt(counts_info$fragments_per_sample))
    } else{
      sampled_expected <-  outer(expectation[background_peaks[peak_set,]], counts_info$fragments_per_sample)
    }
    sampled_deviation = sampled - sampled_expected
  } else{
    tf_vec <- sparseMatrix(j = peak_set, i = rep(1,tf_count), x = 1,
                           dims = c(1, counts_info$npeak))
    
    observed <- as.vector(tf_vec %*% counts_mat)
    if (norm){
      expected <- as.vector(tf_vec %*% (expectation/sqrt(expectation)) %*% (counts_info$fragments_per_sample/sqrt(counts_info$fragments_per_sample)))
    } else{
      expected <- as.vector(tf_vec %*% expectation %*% counts_info$fragments_per_sample)
    }
    observed_deviation = observed - expected
    
    niterations = ncol(background_peaks)
    sample_mat = sparseMatrix(j = as.vector(background_peaks[peak_set,1:niterations]),
                              i = rep(1:niterations, each = tf_count),
                              x=1,
                              dims = c(niterations, counts_info$npeak))
    
    sampled = as.matrix(sample_mat %*% counts_mat);
    if (norm){
      sampled_expected = as.matrix(sample_mat %*% (expectation/sqrt(expectation)) %*% (counts_info$fragments_per_sample/sqrt(counts_info$fragments_per_sample)))
    } else{
      sampled_expected = as.matrix(sample_mat %*% expectation %*% counts_info$fragments_per_sample)
    }
    sampled_deviation = sampled - sampled_expected
    
  }
  
  mean_sampled_deviation <- colMeans(sampled_deviation)
  sd_sampled_deviation <- apply(sampled_deviation, 2, sd)
  
  z <- (observed_deviation - mean_sampled_deviation) / sd_sampled_deviation
  if (length(fail_filter) > 0) z[fail_filter] = NA
  
  v <- sd(z, na.rm = TRUE)
  
  return(v)
}
