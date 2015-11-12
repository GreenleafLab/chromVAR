

compute_variability <- function(motif_indices, counts_mat, bg_param, niterations = 50,
                                metric = c("z-score-binom-2","old","old2","old3","z-score","z-score-2","z-score-3","z-score-binom", "z-score-binom-3"),
                                BPPARAM = BiocParallel::bpparam()){

  metric = match.arg(metric)

  if (class(bg_param) == "deviationBackgroundParameters"){
    results = deviationResultSet(results = BiocParallel::bplapply(motif_indices, compute_deviations0, counts_mat, bg_param, niterations, metric, BPPARAM = BPPARAM), tfs = names(motif_indices))

  } else if (class(bg_param) == "backgroundPeaks"){
    results = deviationResultSet(results = BiocParallel::bplapply(motif_indices, compute_deviations, counts_mat, bg_param, niterations, metric, BPPARAM = BPPARAM), tfs = names(motif_indices))
  } else {
    stop("bg_param argument must be of class deviationBackgroundParameters or backgroundPeaks")
  }
  #results = adjust_p_values(results)

  return(results)
}



compute_deviations0 <- function(peak_set, counts_mat, bg_param, niterations = 50,
                                metric = c("z-score-binom-2","old","old2","old3","z-score","z-score-2","z-score-3","z-score-binom", "z-score-binom-3")){

  metric = match.arg(metric)

  tf_count = length(peak_set)
  tf_vec = sparseMatrix(j = peak_set, i = rep(1,tf_count), x = 1, dims = c(1, length(counts_mat@peaks)))

  observed = as.matrix(tf_vec %*% counts_mat@counts)

  sampled_counts = sampleBackgroundPeaks(object = bg_param, annotation_vector = tf_vec, counts_mat = counts_mat,
                                                   reads = sum(observed), niterations = niterations)

  res <- compute_var_metrics(observed, sampled_counts, counts_mat, metric = metric)

  return(res)
}

compute_deviations <- function(peak_set, counts_mat, bg_param, niterations = 50,
                               metric = c("z-score-binom-2","old","old2","old3","z-score","z-score-2","z-score-3","z-score-binom", "z-score-binom-3")){

  metric = match.arg(metric)

  tf_count = length(peak_set)
  tf_vec = sparseMatrix(j = peak_set, i = rep(1,tf_count), x = 1, dims = c(1, length(counts_mat@peaks)))

  observed = as.matrix(tf_vec %*% counts_mat@counts)

  sampled_counts = sampleBackgroundPeaks(object = bg_param, peak_set = peak_set, counts_mat = counts_mat,
                                                  niterations = niterations)

  res <- compute_var_metrics(observed, sampled_counts, counts_mat, metric = metric)

  return(res)

}



compute_var_metrics <- function(observed, sampled_counts, counts_mat,
                                metric = c("z-score-binom-2","old","old2","old3","z-score","z-score-2","z-score-3","z-score-binom", "z-score-binom-3")){

  metric = match.arg(metric)

  if (metric == 'z-score-binom'){

    expected_prob = sum(observed)/counts_mat@total_fragments
    expected_sampled_prob = colSums(sampled_counts)/counts_mat@total_fragments

    raw_deviation = sapply(1:length(observed), function(x){
      less = ((observed[x]/counts_mat@fragments_per_cell[x])<expected_prob)
      p = pbinom(observed[x], counts_mat@fragments_per_cell[x], prob = expected_prob, lower.tail = less, log.p = TRUE)
      if (less) p = - p
      p
    })

    sampled_deviation = sapply(1:ncol(sampled_counts), function(y) {
      sapply(1:length(observed), function(x){
        less = ((sampled_counts[x,y]/counts_mat@fragments_per_cell[x])<expected_sampled_prob[y])
        p = pbinom(sampled_counts[x,y], counts_mat@fragments_per_cell[x], prob = expected_sampled_prob[y],
                   lower.tail = less, log.p = TRUE)
        if (less) p = - p
        p
      })})

    mean_sampled_deviation = rowMeans(sampled_deviation, na.rm = T)
    sd_sampled_deviation = apply(sampled_deviation, 1, sd, na.rm = T)

    normdev = (raw_deviation - mean_sampled_deviation) / sd_sampled_deviation

    sd_normdev = sd(normdev, na.rm = T)

    bootstrap_indexes = sample(1:length(normdev),length(normdev)*1000,replace=T)
    bootstrap_sds = sapply(1:1000, function(x) sd(normdev[bootstrap_indexes[(1 + (x-1)*length(normdev)):(x*length(normdev))]], na.rm = T))
    sd_error = quantile(bootstrap_sds, c(0.025, 0.975), na.rm=T)

    res = deviationResult(normalized_deviations = normdev, variability = sd_normdev, variability_error = sd_error)

  } else if (metric == 'z-score-binom-2'){

    expected_prob = sum(observed)/counts_mat@total_fragments
    expected_sampled_prob = colSums(sampled_counts)/counts_mat@total_fragments

    expected =  counts_mat@fragments_per_cell * expected_prob
    expected_sampled_counts =  outer(counts_mat@fragments_per_cell, expected_sampled_prob)

    raw_deviation = (observed - expected) / sqrt(expected * (1 - expected_prob))
    sampled_deviation = (sampled_counts - expected_sampled_counts) / sqrt(expected_sampled_counts * (1 - matrix(expected_sampled_prob, nrow = nrow(sampled_counts), ncol = ncol(sampled_counts),byrow = T)))

    mean_sampled_deviation = rowMeans(sampled_deviation)
    sd_sampled_deviation = apply(sampled_deviation, 1, sd)

    normdev = (raw_deviation - mean_sampled_deviation) / sd_sampled_deviation

    sd_normdev = sd(normdev)

    bootstrap_indexes = sample(1:length(normdev),length(normdev)*1000,replace=T)
    bootstrap_sds = sapply(1:1000, function(x) sd(normdev[bootstrap_indexes[(1 + (x-1)*length(normdev)):(x*length(normdev))]]))
    sd_error = quantile(bootstrap_sds, c(0.025, 0.975))

    res = deviationResult(normalized_deviations = as.numeric(normdev), variability = sd_normdev, variability_error = sd_error)

  }else if (metric == 'z-score-binom-3'){

    expected_prob = sum(observed)/counts_mat@total_fragments
    expected_sampled_prob = colSums(sampled_counts)/counts_mat@total_fragments

    expected =  counts_mat@fragments_per_cell * expected_prob
    expected_sampled_counts =  outer(counts_mat@fragments_per_cell, expected_sampled_prob)

    if (min(expected) < 50){

      raw_deviation = sapply(1:length(observed), function(x){
        p = pbinom(observed[x], counts_mat@fragments_per_cell[x], prob = expected_prob, log.p = TRUE)
        qnorm(p, mean = 0, sd = 1, log.p = TRUE)
      })

      sampled_deviation = sapply(1:ncol(sampled_counts), function(y) {
        sapply(1:length(observed), function(x){
          p = pbinom(sampled_counts[x,y], counts_mat@fragments_per_cell[x], prob = expected_sampled_prob[y], log.p = TRUE)
          qnorm(p, mean = 0, sd = 1, log.p = TRUE)
        })})

    } else{

      raw_deviation = (observed - expected) / sqrt(expected * (1 - expected_prob))
      sampled_deviation = (sampled_counts - expected_sampled_counts) / sqrt(expected_sampled_counts * (1 - matrix(expected_sampled_prob, nrow = nrow(sampled_counts), ncol = ncol(sampled_counts),byrow = T)))

    }

    mean_sampled_deviation = rowMeans(sampled_deviation, na.rm = T)
    sd_sampled_deviation = apply(sampled_deviation, 1, sd, na.rm = T)

    normdev = (raw_deviation - mean_sampled_deviation) / sd_sampled_deviation

    sd_normdev = sd(normdev, na.rm = T)

    bootstrap_indexes = sample(1:length(normdev),length(normdev)*1000,replace=T)
    bootstrap_sds = sapply(1:1000, function(x) sd(normdev[bootstrap_indexes[(1 + (x-1)*length(normdev)):(x*length(normdev))]], na.rm = T))
    sd_error = quantile(bootstrap_sds, c(0.025, 0.975), na.rm=T)

    res = deviationResult(normalized_deviations = as.numeric(normdev), variability = sd_normdev, variability_error = sd_error)

  }

  else if (metric == 'z-score'){

    expected =  sum(observed)*counts_mat@fragments_per_cell/counts_mat@total_fragments
    expected_sampled_counts =  outer(counts_mat@fragments_per_cell/counts_mat@total_fragments, colSums(sampled_counts))

    raw_deviation = observed - expected
    mean_sampled_deviation = rowMeans(sampled_counts - expected_sampled_counts)
    sd_sampled_deviation = apply(sampled_counts - expected_sampled_counts, 1, sd)

    normdev = (raw_deviation - mean_sampled_deviation) / sd_sampled_deviation

    sd_normdev = sd(normdev)

    bootstrap_indexes = sample(1:length(normdev),length(normdev)*1000,replace=T)
    bootstrap_sds = sapply(1:1000, function(x) sd(normdev[bootstrap_indexes[(1 + (x-1)*length(normdev)):(x*length(normdev))]]))
    sd_error = quantile(bootstrap_sds, c(0.025, 0.975))

    res = deviationResult(normalized_deviations = as.numeric(normdev), variability = sd_normdev, variability_error = sd_error)

  } else if (metric == 'z-score-2'){

    expected =  sum(observed)*counts_mat@fragments_per_cell/counts_mat@total_fragments
    expected_sampled_counts =  outer(counts_mat@fragments_per_cell/counts_mat@total_fragments, colSums(sampled_counts))

    raw_deviation = (observed - expected) / expected
    sampled_deviation = (sampled_counts - expected_sampled_counts) / expected_sampled_counts
    mean_sampled_deviation = rowMeans(sampled_deviation)
    sd_sampled_deviation = apply(sampled_deviation, 1, sd)

    normdev = (raw_deviation - mean_sampled_deviation) / sd_sampled_deviation

    sd_normdev = sd(normdev)

    bootstrap_indexes = sample(1:length(normdev),length(normdev)*1000,replace=T)
    bootstrap_sds = sapply(1:1000, function(x) sd(normdev[bootstrap_indexes[(1 + (x-1)*length(normdev)):(x*length(normdev))]]))
    sd_error = quantile(bootstrap_sds, c(0.025, 0.975))

    res = deviationResult(normalized_deviations = as.numeric(normdev), variability = sd_normdev, variability_error = sd_error)

  } else if (metric == 'z-score-3'){

    expected =  sum(observed)*counts_mat@fragments_per_cell/counts_mat@total_fragments
    expected_sampled_counts =  outer(counts_mat@fragments_per_cell/counts_mat@total_fragments, colSums(sampled_counts))

    pseudo = 50
    raw_deviation = (observed - expected) / (expected + pseudo)
    sampled_deviation = (sampled_counts - expected_sampled_counts) / (expected_sampled_counts + pseudo)
    mean_sampled_deviation = rowMeans(sampled_deviation)
    sd_sampled_deviation = apply(sampled_deviation, 1, sd)

    normdev = (raw_deviation - mean_sampled_deviation) / sd_sampled_deviation

    sd_normdev = sd(normdev)

    bootstrap_indexes = sample(1:length(normdev),length(normdev)*1000,replace=T)
    bootstrap_sds = sapply(1:1000, function(x) sd(normdev[bootstrap_indexes[(1 + (x-1)*length(normdev)):(x*length(normdev))]]))
    sd_error = quantile(bootstrap_sds, c(0.025, 0.975))

    res = deviationResult(normalized_deviations = as.numeric(normdev), variability = sd_normdev, variability_error = sd_error)

  } else if (metric == 'old3'){

    expected =  sum(observed)*counts_mat@fragments_per_cell/counts_mat@total_fragments
    expected_sampled_counts =  outer(counts_mat@fragments_per_cell/counts_mat@total_fragments, colSums(sampled_counts))

    raw_deviation = observed - expected
    mean_sampled_deviation = rowMeans(sampled_counts - expected_sampled_counts)
    rms_sampled_deviation = apply(sampled_counts - expected_sampled_counts, 1, rms)

    normdev = (raw_deviation - mean_sampled_deviation) / rms_sampled_deviation

    sd_normdev = sd(normdev)

    bootstrap_indexes = sample(1:length(normdev),length(normdev)*1000,replace=T)
    bootstrap_sds = sapply(1:1000, function(x) sd(normdev[bootstrap_indexes[(1 + (x-1)*length(normdev)):(x*length(normdev))]]))
    sd_error = quantile(bootstrap_sds, c(0.025, 0.975))

    res = deviationResult(normalized_deviations = as.numeric(normdev), variability = sd_normdev, variability_error = sd_error)

  } else if (metric == 'old'){

    expected =  sum(observed)*counts_mat@fragments_per_cell/counts_mat@total_fragments
    expected_sampled_counts =  outer(counts_mat@fragments_per_cell/counts_mat@total_fragments, colSums(sampled_counts))

    raw_deviation = observed - expected
    sampled_deviation = sampled_counts-expected_sampled_counts
    rms_sampled_deviation = apply(sampled_counts - expected_sampled_counts, 1, rms)

    normdev = raw_deviation / rms_sampled_deviation

    normvar_func <- function(raw_deviation, sampled_deviation){
      sqrt(sum(raw_deviation**2)/sum(rowMeans(sampled_deviation**2)))
    }

    normvar = normvar_func(raw_deviation, sampled_deviation)

    bootstrap_indexes = sample(1:length(normdev),length(normdev)*1000,replace=T)
    bootstrap_normvars = sapply(1:1000, function(x)
      normvar_func(raw_deviation[bootstrap_indexes[(1 + (x-1)*length(normdev)):(x*length(normdev))]],sampled_deviation[bootstrap_indexes[(1 + (x-1)*length(normdev)):(x*length(normdev))],]))
    normvar_error = quantile(bootstrap_normvars, c(0.025, 0.975))

    res = deviationResult(normalized_deviations = as.numeric(normdev), variability = normvar, variability_error = normvar_error)

  } else if (metric == 'old2'){

    expected =  sum(observed)*counts_mat@fragments_per_cell/counts_mat@total_fragments
    expected_sampled_counts =  outer(counts_mat@fragments_per_cell/counts_mat@total_fragments, colSums(sampled_counts))

    raw_deviation = observed - expected
    rms_sampled_deviation = apply(sampled_counts - expected_sampled_counts, 1, rms)

    normdev = raw_deviation / rms_sampled_deviation

    sd_normdev = sd(normdev)

    bootstrap_indexes = sample(1:length(normdev),length(normdev)*1000,replace=T)
    bootstrap_sds = sapply(1:1000, function(x) sd(normdev[bootstrap_indexes[(1 + (x-1)*length(normdev)):(x*length(normdev))]]))
    sd_error = quantile(bootstrap_sds, c(0.025, 0.975))

    res = deviationResult(normalized_deviations = as.numeric(normdev), variability = sd_normdev, variability_error = sd_error)

  }
  return(res)
}
