

pwm_to_prob_helper <- function(pwm){
  type <- pwm_type(pwm)
  mat <- TFBSTools::as.matrix(pwm)
  if (type == "frequency"){
    out <- mat
  } else if (type == "log2"){
    out <- 2**mat * matrix(TFBSTools::bg(pwm), nrow = nrow(mat), ncol= ncol(mat), byrow = FALSE)
  } else if (type == "log"){
    out <- exp(mat)* matrix(TFBSTools::bg(pwm), nrow = nrow(mat), ncol= ncol(mat), byrow = FALSE)
  }
  return(out)
}



pwm_to_prob <- function(pwms){
  if (inherits(pwms, "PWMatrix")){
    out <- list(pwm_to_prob_helper(pwms))
  } else if (inherits(pwms, "PWMatrixList")){
    out <- lapply(pwms, pwm_to_prob_helper)
  } else if (inherits(pwms,"matrix")){
    stopifnot(all.equal(colSums(pwms),rep(1,ncol(pwms))))
    out <- list(pwms)
  } else if (inherits(pwms, "list")){
    stopifnot(all_true(sapply(pwms, function(x) all.equal(colSums(x), rep(1,ncol(x))))))
    out <- pwms
  } else{
    stop("incorrect input format, must be PWMatrix,PWMatrixList, matrix, list of matrices")
  }
  return(out)
}



#' pwm_distance
#'
#' @param x Feature 1
#' @param y Feature 2
#' @param min_overlap 
#'
#' @return A distance metric between the PWM
#' @export
pwm_distance <- function(x, y = NULL, min_overlap = 5){
  if (is.null(y)){
    mats <- pwm_to_prob(x)
    compute_pwm_dist(mats, min_overlap)
  } else{
    mats1 <- pwm_to_prob(x)
    mats2 <- pwm_to_prob(y)
    compute_pwm_dist2(mats1, mats2, min_overlap)
  }
}
