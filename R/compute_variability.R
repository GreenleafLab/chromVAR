
quantile_helper <- function(values, quantiles, na.rm) {
  if (na.rm) {
    out <- quantile(values, quantiles, na.rm = TRUE)
  } else {
    if (!all_false(is.na(values))) {
      out <- rep(NA, length(quantiles))
    } else {
      out <- quantile(values, quantiles)
    }
  }
  return(out)
}

#'computeVariability
#'
#'function to compute overall variability of motif sets across samples
#'@param object output from \code{\link{computeDeviations}}
#'@param bootstrap_error compute bootstrap confidence interval
#'@param bootstrap_samples number of bootstrap samples to take
#'@param bootstrap_quantiles quantiles for bootstrap 
#'@param na.rm remove NAs? default is true
#'@return data.frame with columns for name, variability, bootstrap lower bound,
#'bootstrap upper bound, raw p value, adjust p value.
#'@export
#'@examples
#' # Load very small example results from computeDeviations
#' data(mini_dev, package = "chromVAR")
#' variability <- computeVariability(mini_dev)                        
computeVariability <- function(object, 
                                bootstrap_error = TRUE, 
                                bootstrap_samples = 1000, 
                                bootstrap_quantiles = c(0.025, 0.975), 
                                na.rm = TRUE) {
  
  stopifnot(inherits(object, "chromVARDeviations") || 
              canCoerce(object, "chromVARDeviations"))
  
  sd_deviations <- row_sds(assays(object)$z, na.rm)
  
  p_sd <- pchisq((ncol(object) - 1) * (sd_deviations^2), 
                 df = (ncol(object) - 1), 
                 lower.tail = FALSE)
  
  p_adj <- p.adjust(p = p_sd, method = "BH")
  
  name_vals <- if (!is.null(rowData(object)$name)) 
    rowData(object)$name else rownames(object)
  
  if (bootstrap_error) {
    # check some arguments
    stopifnot(bootstrap_samples > 0 && all_whole(bootstrap_samples))
    stopifnot(all(bootstrap_quantiles > 0))
    stopifnot(all(bootstrap_quantiles < 1))
    stopifnot(bootstrap_quantiles[2] > bootstrap_quantiles[1])
    
    boot_sd <- do.call(rbind, bplapply(seq_len(bootstrap_samples), 
                                       function(x) 
                                         row_sds_perm(assays(object)$z, 
                                                                na.rm)))
    sd_error <- apply(boot_sd, 2, quantile_helper, 
                      quantiles = bootstrap_quantiles, 
                      na.rm = na.rm)
    
    out <- data.frame(name = name_vals, 
                      variability = sd_deviations, 
                      bootstrap_lower_bound = sd_error[1, ], 
                      bootstrap_upper_bound = sd_error[2, ],
                      p_value = p_sd, 
                      p_value_adj = p_adj, 
                      row.names = rownames(object))
  } else {
    out <- data.frame(name = name_vals,
                      variability = sd_deviations, 
                      p_value = p_sd, 
                      p_value_adj = p_adj, 
                      row.names = rownames(object))
  }
  
  return(out)
}



# Helper function --------------------------------------------------------------

compute_variability_single <- function(peak_set, 
                                       counts_mat, 
                                       background_peaks, 
                                       expectation) {
  
  tmp <- compute_deviations_single(peak_set, counts_mat, background_peaks,
                                   expectation)
  
  v <- sd(tmp$z, na.rm = TRUE)
  
  return(v)
}
