pwm_type <- function(pwm) {
  # Determine whether un-logged, natural log, or log2
  if (isTRUE(all.equal(colSums(as.matrix(pwm)), 
                       rep(1, length(pwm)),
                       check.attributes = FALSE))) {
    return("frequency")
  } else if (isTRUE(all.equal(colSums(2^(as.matrix(pwm)) *
                                      matrix(bg(pwm),
                                             byrow = FALSE,
                                             ncol = length(pwm),
                                             nrow = 4)),
                              rep(1, length(pwm)), tolerance = 10^-5,
                              check.attributes = FALSE))) {
    return("log2")
  } else if (isTRUE(all.equal(colSums(exp(as.matrix(pwm)) *
                                      matrix(bg(pwm),
                                             byrow = FALSE,
                                             ncol = length(pwm),
                                             nrow = 4)),
                              rep(1, length(pwm)), tolerance = 10^-5,
                              check.attributes = FALSE))) {
    return("log")
  } else {
    stop("Can't determine format of PWM -- should be numeric ",
         "frequency summing to 1 or log or log2 odds ratio")
  }
}

pwm_to_prob_helper <- function(pwm) {
  type <- pwm_type(pwm)
  mat <- TFBSTools::as.matrix(pwm)
  if (type == "frequency") {
    out <- mat
  } else if (type == "log2") {
    out <- 2^mat * matrix(TFBSTools::bg(pwm), nrow = nrow(mat), 
                          ncol = ncol(mat), 
      byrow = FALSE)
  } else if (type == "log") {
    out <- exp(mat) * matrix(TFBSTools::bg(pwm), nrow = nrow(mat), 
                             ncol = ncol(mat), 
      byrow = FALSE)
  }
  return(out)
}

pfm_to_prob_helper <- function(pfm){
  mat <- TFBSTools::as.matrix(pfm)
  mat / matrix(colSums(mat), byrow = TRUE, nrow = nrow(mat), ncol = ncol(mat))
}

pwm_to_prob <- function(pwms) {
  if (inherits(pwms, "PWMatrix")) {
    out <- list(pwm_to_prob_helper(pwms))
  } else if (inherits(pwms, "PWMatrixList")) {
    out <- lapply(pwms, pwm_to_prob_helper)
  } else if (inherits(pwms, "matrix")) {
    stopifnot(all.equal(colSums(pwms), rep(1, ncol(pwms))))
    out <- list(pwms)
  } else if (inherits(pwms, "PFMatrix")) {
    out <- list(pfm_to_prob_helper(pwms))
  } else if (inherits(pwms, "PFMatrixList")) {
    out <- lapply(pwms, pfm_to_prob_helper)
  }else if (inherits(pwms, "list")) {
    stopifnot(all(vapply(pwms, function(x){
      all(abs(colSums(x) - 1) < 10e-5)},
      TRUE)))
    out <- pwms
  } else {
    stop("incorrect input format, must be PWMatrix,PWMatrixList, matrix, 
         list of matrices")
  }
  return(out)
}



#' pwmDistance
#'
#' computes distance between every pwm in a list or between pwms in one list
#' with pwms in another
#' @param x list of pwms or pfms, see Details
#' @param y list of pwms or pfms, see Details
#' @param min_overlap minimum number of basepairs overlapping between motifs
#' @details The format of x and y should be a 
#' \code{\link[TFBSTools]{PWMatrixList}} 
#' or \code{\link[TFBSTools]{PFMatrixList}} or a list of matrices with rows 
#' corresponding to "A","C","G","T" and columns summing to 1. 
#' @return a list with three matrices- 'dist' has the distance between each
#'  pair of motifs, 'strand' has the strand of the motif for the match, and 
#'  'offset' has the offset between the motifs. 
#' @export
#' @examples 
#' 
#' motifs <- getJasparMotifs()
#' library(TFBSTools)
#' pwm_dists <- pwmDistance(toPWM(motifs[[1]]), toPWM(motifs[[2]]))
pwmDistance <- function(x, y = NULL, min_overlap = 5) {
  if (is.null(y)) {
    mats <- pwm_to_prob(x)
    compute_pwm_dist(mats, min_overlap)
  } else {
    mats1 <- pwm_to_prob(x)
    mats2 <- pwm_to_prob(y)
    compute_pwm_dist2(mats1, mats2, min_overlap)
  }
}
