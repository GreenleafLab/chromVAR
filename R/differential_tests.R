#' differential_deviations
#'
#' Function to see whether deviations differ between groups
#' @param object chromVARDeviations object
#' @param groups either vector of groups or name of column in colData of object with grouop information
#' @param alternative only used if there are two groups -- two.sided or one sided test
#' @param parametric use parametric test. alternatively will use kruskal wallace
#' @return data.frame with p value and adjusted p value
#' @export
differential_deviations <- function(object, groups,
                                       alternative = c("two.sided", "less",
                                                       "greater"), parametric = TRUE) {
  stopifnot(is(object,"chromVARDeviations"))
  if (length(groups) == 1 && groups %in% colnames(colData(object))){
    groups <- colData(object)[groups]
  } else if (length(groups) != ncol(object)){
    stop("invalid groups input, must be vector of lench ncol(object) or column name from colData(object)")
  }

  inputs <- lapply(split(1:ncol(object),groups), function(x) deviations(object)[,x])
  alternative <- match.arg(alternative)
  if (parametric) {
    if (length(inputs) == 2) {
      # t-test
      p_val <- sapply(1:nrow(inputs[[1]]),
                      function(x) t_helper(inputs[[1]][x,],
                                           inputs[[2]][x, ], 
                                           alternative))
    } else {
      # anova
      reshaped <- lapply(1:nrow(inputs[[1]]),
                         function(x) lapply(inputs, function(y) y[x,]))
      p_val <- sapply(reshaped, function(x) do.call(anova_helper, x))
    }
  } else {
    if (length(inputs) == 2) {
      # wilcoxon
      p_val <- sapply(1:nrow(inputs[[1]]),
                      function(x) wilcoxon_helper(inputs[[1]][x, ],
                                                  inputs[[2]][x, ],
                                                  alternative))
    } else {
      # kruskal-wallis
      reshaped <- lapply(1:nrow(inputs[[1]]),
                         function(x) lapply(inputs, function(y) y[x,
                                                                  ]))
      p_val <- sapply(reshaped, function(x) do.call(kw_helper, x))
    }
  }
  p_adj <- p.adjust(p_val, method = "BH")
  return(data.frame(p_value = p_val, p_value_adjusted = p_adj))
}

kw_helper <- function(...) {
  inputs <- list(...)
  vals <- do.call(c, inputs)
  group <- factor(do.call(c, lapply(1:length(inputs),
                                    function(x) rep(x, length(inputs[[x]])))))
  res <- kruskal.test(vals ~ group)
  return(res$p.value)
}

wilcoxon_helper <- function(x, y, alternative) {
  return(wilcox.test(x, y, alternative = alternative, paired = FALSE)$p.value)
}


t_helper <- function(x, y, alternative) {
  return(t.test(x, y, alternative = alternative, paired = FALSE,
                var.equal = FALSE)$p.value)
}

anova_helper <- function(...) {
  inputs <- list(...)
  vals <- do.call(c, inputs)
  group <- factor(do.call(c, lapply(1:length(inputs),
                                    function(x) rep(x, length(inputs[[x]])))))
  anova_res <- oneway.test(vals ~ group, var.equal = FALSE)
  return(anova_res$p.value)
}


#' differential_variability
#'
#' Function to determine whether groups differ in variability
#' @param object chromVARDeviations object
#' @param groups either vector of groups or name of column in colData of object with grouop information
#' @param parametric use parametric test. alternatively will use kruskal wallace
#' @return data.frame with p value and adjusted p value
#' @export
differential_variability <- function(object, groups, parametric = TRUE) {
  stopifnot(is(object,"chromVARDeviations"))
  if (length(groups) == 1 && groups %in% colnames(colData(object))){
    groups <- colData(object)[groups]
  } else if (length(groups) != ncol(object)){
    stop("invalid groups input, must be vector of lench ncol(object) or column name from colData(object)")
  }
  # Brown-Forsythe test
  inputs <- lapply(split(1:ncol(object),groups),
                   function(x) deviation_scores(object)[,x])
  reshaped <- lapply(1:nrow(inputs[[1]]),
                     function(x) lapply(inputs, function(y) y[x, ]))
  if (parametric) {
    p_val <- sapply(reshaped, function(x) do.call(bf_var_test, x))
  } else {
    p_val <- sapply(reshaped, function(x) do.call(bf_kw_var_test, x))
  }
  p_adj <- p.adjust(p_val, method = "BH")
  return(data.frame(p_value = p_val, p_value_adjusted = p_adj))
}

bf_var_test <- function(...) {
  inputs <- list(...)
  medians <- sapply(inputs, median, na.rm = TRUE)
  median_diff <- do.call(c, lapply(1:length(inputs),
                                   function(x) abs(inputs[[x]] -
                                                     medians[x])))
  group <- factor(do.call(c, lapply(1:length(inputs),
                                    function(x) rep(x, length(inputs[[x]])))))
  return(anova(lm(median_diff ~ group))[1, 5])
}

bf_kw_var_test <- function(...) {
  inputs <- list(...)
  medians <- sapply(inputs, median, na.rm = TRUE)
  median_diff <- do.call(c, lapply(1:length(inputs),
                                   function(x) abs(inputs[[x]] -
                                                     medians[x])))
  group <- factor(do.call(c, lapply(1:length(inputs),
                                    function(x) rep(x, length(inputs[[x]])))))
  return(res <- kruskal.test(median_diff ~ group)$p.value)
}

