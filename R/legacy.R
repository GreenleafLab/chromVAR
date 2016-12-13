
compute_deviations_legacy <- function(object,
                                      annotations = NULL,
                                      background_peaks = NULL,
                                      expectation = NULL){


  stopifnot(inherits(object, "SummarizedExperiment"))

  if (is.null(assays(object)$counts)) stop("No counts slot")

  if (inherits(assays(object)$counts,"matrix")){
    assays(object)$counts = Matrix(counts_mat)
  }
  stopifnot(inherits(assays(object)$counts,"Matrix"))

  if (min(get_fragments_per_peak(object)) <= 0) stop("All peaks must have at least one fragment in one sample")

  if (is.null(background_peaks)){
    background_peaks <- get_background_peaks(object)
  }

  stopifnot(nrow(object) == nrow(background_peaks))

  if (inherits(annotations, "SummarizedExperiment")){
    peak_indices <- convert_to_ix_list(assays(annotations)$match)
  } else {
    stop("peak_indices must be given as SummarizedExperiment with a 'match' slot")
  }

  if (is.null(expectation)){
    expectation <- compute_expectations(object)
  } else{
    stopifnot(length(expectation) == nrow(object))
  }

  # check that indices fall within appropriate bounds
  tmp <- unlist(peak_indices, use.names =F)
  if (is.null(tmp) ||
      !(all.equal(tmp, as.integer(tmp))) ||
      max(tmp) > nrow(object) ||
      min(tmp) < 1){
    stop("peak_indices are not valid")
  }

  if(is.null(names(peak_indices))){
    names(peak_indices) = 1:length(peak_indices)
  }

  sample_names <- colnames(object)

  results <- BiocParallel::bplapply(peak_indices,
                                    compute_deviations_single_legacy,
                                    counts_mat = assays(object)$counts,
                                    background_peaks = background_peaks,
                                    expectation = expectation)


  dev <- t(vapply(results, function(x) x[["deviations"]], rep(0,ncol(object))))

  colnames(dev) = sample_names

  out <- SummarizedExperiment(assays = list(deviations = dev),
                              colData = colData(object),
                              rowData = cbind(colData(annotations),
                                              DataFrame(variability = sapply(results, function(x) x[["variability"]]))))

  return(out)
}



compute_deviations_single_legacy <- function(peak_set,
                                      counts_mat,
                                      background_peaks,
                                      expectation,
                                      intermediate_results = FALSE){


  tf_count <- length(peak_set)
  fragments_per_sample = colSums(counts_mat)

  if (tf_count == 1){
    observed <- as.vector(counts_mat[peak_set,])
    names(observed) <- colnames(counts_mat)
    sampled_counts <-  t(as.matrix(counts_mat[background_peaks[peak_set,],]))
  }
  else {
    tf_vec <- sparseMatrix(j = peak_set, i = rep(1,tf_count), x = 1,
                           dims = c(1, nrow(counts_mat)))

    observed <- as.vector(tf_vec %*% counts_mat)
    names(observed) <- colnames(counts_mat)

    sampled_counts <- sample_background_peaks_legacy(counts_mat = counts_mat,
                                              background_peaks = background_peaks,
                                              peak_set = peak_set)
  }


  expected <-  fragments_per_sample * sum(expectation[peak_set])
  expected_sampled_counts <- outer(apply(background_peaks, 2, function(x) sum(expectation[x[peak_set]])),
                                   fragments_per_sample)


  observed_deviation <- observed - expected
  sampled_deviation <- sampled_counts - expected_sampled_counts
  rms_sampled_deviation <- apply(sampled_deviation, 2, rms)
  normdev <- observed_deviation / rms_sampled_deviation
  normvar <- sqrt(sum(observed_deviation**2)/mean(rowSums(sampled_deviation**2)))
  if (intermediate_results){
      out = list(deviations = normdev, variability = normvar, observed = observed,
                 sampled = sampled, expected = expected, expected_sampled = expected_sampled_counts)
  } else{
    out = list(deviations = normdev, variability = normvar)
  }
  return(out)
}


sample_background_peaks_legacy <- function(counts_mat, peak_set, background_peaks){

  niterations = ncol(background_peaks)
  sample_mat = sparseMatrix(j = as.vector(background_peaks[peak_set,1:niterations]),
                            i = rep(1:niterations, each = length(peak_set)),
                            x=1,
                            dims = c(ncol(background_peaks), nrow(counts_mat)))

  sampled_counts =  as.matrix(sample_mat %*% counts_mat)

  return(sampled_counts)
}



